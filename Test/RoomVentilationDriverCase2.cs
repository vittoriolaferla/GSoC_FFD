using System;
using System.Collections.Generic;
using System.Globalization; // For culture settings
using System.IO;          // For file and directory operations
using System.Numerics;    // For Vector3
using System.Text.RegularExpressions;
using System.Threading;
using MathNet.Numerics.Optimization;

namespace FastFluidSolver
{
    /// <summary>
    /// Driver for the room ventilation simulation.
    /// </summary>
    class RoomVentilationDriverCase2
    {
        // Class to hold vent parameters read from the CSV file.
        public class VentParameter
        {
            public string vent_type;
            public string surf_id;
            public double XB1;
            public double XB2;
            public double XB3;
            public double XB4;
            public double XB5;
            public double XB6;
            public double VEL;
            public string VEL_T; // For louvered and swirl vents, this field should contain two comma-separated values.
        }

        // Reads the CSV file (using a regex to handle quoted fields) and returns a list of vent parameters.
        static List<VentParameter> LoadVentParameters(string filename)
        {
            var ventList = new List<VentParameter>();

            if (!File.Exists(filename))
            {
                Console.WriteLine($"Vent parameter file {filename} not found.");
                return ventList;
            }

            var lines = File.ReadAllLines(filename);
            // Use Regex.Split to split on commas that are not inside quotes.
            for (int i = 1; i < lines.Length; i++)
            {
                string line = lines[i].Trim();
                if (string.IsNullOrEmpty(line))
                    continue;
                // This regex splits on commas only if they are followed by an even number of quotes.
                var tokens = Regex.Split(line, ",(?=(?:[^\"]*\"[^\"]*\")*[^\"]*$)");
                if (tokens.Length < 9)
                    continue;

                var vent = new VentParameter();
                vent.vent_type = tokens[0].Trim();
                vent.surf_id = tokens[1].Trim();
                vent.XB1 = double.Parse(tokens[2], CultureInfo.InvariantCulture);
                vent.XB2 = double.Parse(tokens[3], CultureInfo.InvariantCulture);
                vent.XB3 = double.Parse(tokens[4], CultureInfo.InvariantCulture);
                vent.XB4 = double.Parse(tokens[5], CultureInfo.InvariantCulture);
                vent.XB5 = double.Parse(tokens[6], CultureInfo.InvariantCulture);
                vent.XB6 = double.Parse(tokens[7], CultureInfo.InvariantCulture);
                vent.VEL = double.Parse(tokens[8], CultureInfo.InvariantCulture);
                if (tokens.Length >= 10)
                    vent.VEL_T = tokens[9].Trim(); // This field will be something like "-0.0707106781186548,0.0"
                ventList.Add(vent);
            }
            return ventList;
        }


        static void RunSimulationForCSV(string ventParamFile, string parentOutputDir)
        {
            Console.WriteLine("Running simulation for CSV: " + ventParamFile);

            // Simulation parameters.
            int Nx = 47; // Number of cells in x (excluding ghost cells)
            int Ny = 47; // Number of cells in y (excluding ghost cells)
            int Nz = 47; // Number of cells in z (excluding ghost cells)
            double length_x = 4.0, length_y = 5.0, length_z = 3.0;
            double dt = 1, nu = 1e-2;
            double tf = 5, t = 0;

            // Set initial conditions.
            double[,,] u0 = new double[Nx + 1, Ny + 2, Nz + 2];
            double[,,] v0 = new double[Nx + 2, Ny + 1, Nz + 2];
            double[,,] w0 = new double[Nx + 2, Ny + 2, Nz + 1];

            double[,,] f_x = new double[Nx + 1, Ny + 2, Nz + 2];
            double[,,] f_y = new double[Nx + 2, Ny + 1, Nz + 2];
            double[,,] f_z = new double[Nx + 2, Ny + 2, Nz + 1];

            // Solver parameters.
            FluidSolver.solver_struct solver_prams = new FluidSolver.solver_struct
            {
                tol = 1e-4,
                min_iter = 0,
                max_iter = 30,
                verbose = true,
                backtrace_order = 2
            };

            // Create the domain including ghost cells.
            RoomVentilationDomain omega = new RoomVentilationDomain(
                Nx + 2, Ny + 2, Nz + 2, length_x, length_y, length_z);

            // Load vent parameters from the CSV file.
            var vents = LoadVentParameters(ventParamFile);

            // Create a subdirectory for this CSV file's simulation.
            string csvName = Path.GetFileNameWithoutExtension(ventParamFile);
            string outputDirectory = Path.Combine(parentOutputDir, csvName);
            Directory.CreateDirectory(outputDirectory);

            // Create subdirectories for final and intermediate output.
            string finalDirectory = Path.Combine(outputDirectory, "final");
            string intermediateDirectory = Path.Combine(outputDirectory, "intermediate");
            Directory.CreateDirectory(finalDirectory);
            Directory.CreateDirectory(intermediateDirectory);

            // Loop over each vent parameter and instantiate the vent.
            foreach (var vent in vents)
            {
                // Create start and end positions based on the XB parameters.
                Vector3 start = new Vector3((float)vent.XB1, (float)vent.XB3, (float)vent.XB5);
                Vector3 end = new Vector3((float)vent.XB2, (float)vent.XB4, (float)vent.XB6);

                // Determine the total velocity magnitude based on surf_id.
                float totalVelocity;
                if (vent.surf_id.IndexOf("supply", StringComparison.OrdinalIgnoreCase) >= 0)
                    totalVelocity = (float)Math.Abs(vent.VEL);
                else if (vent.surf_id.IndexOf("exhaust", StringComparison.OrdinalIgnoreCase) >= 0)
                    totalVelocity = -(float)Math.Abs(vent.VEL);
                else
                    totalVelocity = (float)vent.VEL;

                // If the vent is an exhaust, always use the linear (standard) inlet.
                if (vent.surf_id.IndexOf("exhaust", StringComparison.OrdinalIgnoreCase) >= 0)
                {
                    omega.AddStandardInlet(start, end, totalVelocity);
                }
                // For supply vents: if the vent type is Louvered or Swirl, treat it as a louvered vent.
                else if (vent.vent_type.StartsWith("Louvered Square", StringComparison.OrdinalIgnoreCase) ||
                         vent.vent_type.StartsWith("Swirl", StringComparison.OrdinalIgnoreCase))
                {
                    // For louvered and swirl vents, the VEL_T field should contain two comma-separated values.
                    string vtStr = vent.VEL_T.Replace("\"", "").Trim();
                    var vtTokens = vtStr.Split(new char[] { ',' }, StringSplitOptions.RemoveEmptyEntries);
                    if (vtTokens.Length != 2)
                    {
                        Console.WriteLine("Warning: Invalid VEL_T parameter for vent: " + vtStr);
                        return;  // Abort this simulation and go to the next CSV.
                    }
                    float vT1 = float.Parse(vtTokens[0], CultureInfo.InvariantCulture);
                    float vT2 = float.Parse(vtTokens[1], CultureInfo.InvariantCulture);
                    omega.AddLouveredInlet(start, end, vT1, vT2, totalVelocity);
                }
                else if (vent.vent_type.StartsWith("Linear", StringComparison.OrdinalIgnoreCase))
                {
                    omega.AddStandardInlet(start, end, totalVelocity);
                }
                else
                {
                    Console.WriteLine($"Unknown vent type: {vent.vent_type}. Aborting simulation for CSV: {ventParamFile}.");
                    return;  // Abort this simulation and go to the next CSV.
                }
            }

            // Initialize the fluid solver and post processor.
            FluidSolver ffd = new FluidSolver(omega, dt, nu, u0, v0, w0, solver_prams);
            PostProcessor pp = new PostProcessor(ffd, omega);

            int tstep = 0;
            // Export initial data and geometry to the final directory.
            pp.export_data_vtk(Path.Combine(finalDirectory, $"roomVentilation_{tstep}.vtk"), 0, false);
            pp.export_geometry_vtk(Path.Combine(finalDirectory, "roomVentilation_geometry.vtk"), 0);

            // Time-stepping loop.
            while (t < tf)
            {
                t += dt;
                tstep++;

                Console.WriteLine("Time t = {0}", t);
                ffd.time_step(f_x, f_y, f_z);

                // 1) Export final velocity & pressure to the final directory:
                pp.export_data_vtk(Path.Combine(finalDirectory,
                   $"roomVentilation_{tstep}.vtk"), t, false);

                // 2) Export the intermediate velocity to the intermediate directory:
                pp.export_intermediate_vtk(Path.Combine(intermediateDirectory,
                   $"roomVentilation_intermediate_{tstep}.vtk"), t, false);
            }

            Console.WriteLine("Simulation for CSV " + ventParamFile + " completed.");
        }

        static void Main()
        {
            // Set the culture.
            Thread.CurrentThread.CurrentCulture = CultureInfo.InvariantCulture;
            Thread.CurrentThread.CurrentUICulture = CultureInfo.InvariantCulture;

            // Define the base directory as "../../../../../" relative to the current working directory.
            string baseDirectorySimulation = Path.GetFullPath(Path.Combine(Directory.GetCurrentDirectory(), "../../../../../../"));

             // Define the base directory as "../../../../../" relative to the current working directory.
            string baseDirectoryCSV = Path.GetFullPath(Path.Combine(Directory.GetCurrentDirectory(), "../../../../../"));

            // Create a new parent directory to store all VTK files from this run under the base directory.
            // The folder name includes a timestamp so that each run is stored separately.
            string parentOutputDir = Path.Combine(
                baseDirectorySimulation,
                "simulation_output",
                DateTime.Now.ToString("yyyyMMdd_HHmmss"));
            Directory.CreateDirectory(parentOutputDir);
            Console.WriteLine("Output directory for simulation run: " + parentOutputDir);

            // Specify the folder containing CSV files, assumed to be at "../../../../../csv_files/"
            string csvFolder = Path.Combine(baseDirectoryCSV, "csv_files");
            string[] csvFiles = Directory.GetFiles(csvFolder, "*.csv");

            // Loop over each CSV file and run the simulation.
            foreach (var csvFile in csvFiles)
            {
                RunSimulationForCSV(csvFile, parentOutputDir);
            }

            Console.WriteLine("All simulations completed.");
        }
    }
}
