using System;
using System.Collections.Generic;
using System.Numerics;

namespace FastFluidSolver
{
    /// <summary>
    /// Domain for room ventilation simulations with variable physical sizes and 2D inlets fully attached to walls.
    /// </summary>
    public class RoomVentilationDomain : Domain
    {
        private List<Inlet> inlets;
        private readonly VelocityCalculator velocityCalculator;
        private readonly PhysicalToCellMapper cellMapper;


        /// <summary>
        /// Constructor with variable domain lengths
        /// </summary>
        /// <param name="Nx">Number of cells in x direction (including ghost cells)</param>
        /// <param name="Ny">Number of cells in y direction (including ghost cells)</param>
        /// <param name="Nz">Number of cells in z direction (including ghost cells)</param>
        /// <param name="length_x">Physical length of the domain in x direction</param>
        /// <param name="length_y">Physical length of the domain in y direction</param>
        /// <param name="length_z">Physical length of the domain in z direction</param>
        public RoomVentilationDomain(int Nx, int Ny, int Nz, double length_x, double length_y, double length_z)
        {
            this.Nx = Nx;
            this.Ny = Ny;
            this.Nz = Nz;

            this.length_x = length_x;
            this.length_y = length_y;
            this.length_z = length_z;

            // Calculate cell sizes
            hx = length_x / (Nx - 2);
            hy = length_y / (Ny - 2);
            hz = length_z / (Nz - 2);

            boundary_cells = new int[Nx, Ny, Nz];
            obstacle_cells = new int[Nx, Ny, Nz];
            boundary_normal_x = new int[Nx, Ny, Nz];
            boundary_normal_y = new int[Nx, Ny, Nz];
            boundary_normal_z = new int[Nx, Ny, Nz];

            boundary_u = new double[Nx - 1, Ny, Nz];
            boundary_v = new double[Nx, Ny - 1, Nz];
            boundary_w = new double[Nx, Ny, Nz - 1];

            outflow_boundary_x = new int[Nx, Ny, Nz];
            outflow_boundary_y = new int[Nx, Ny, Nz];
            outflow_boundary_z = new int[Nx, Ny, Nz];


            // Initialize boundary and obstacle flags
            set_ghost_flags();
            set_boundary_flags();

            inlets = new List<Inlet>();
            velocityCalculator = new VelocityCalculator();

            // The mapper requires Nx, Ny, Nz that include ghost cells in many codes,
            // but you can adapt if you need different indexing.
            cellMapper = new PhysicalToCellMapper(Nx, Ny, Nz, length_x, length_y, length_z);
        }

        /// <summary>
        /// Copy constructor with deep cloning of arrays
        /// </summary>
        /// <param name="old">Existing RoomVentilationDomain to copy from</param>
        public RoomVentilationDomain(RoomVentilationDomain old)
        {
            Nx = old.Nx;
            Ny = old.Ny;
            Nz = old.Nz;

            hx = old.hx;
            hy = old.hy;
            hz = old.hz;

            length_x = old.length_x;
            length_y = old.length_y;
            length_z = old.length_z;

            boundary_cells = old.boundary_cells;
            obstacle_cells = old.obstacle_cells;
            boundary_normal_x = old.boundary_normal_x;
            boundary_normal_y = old.boundary_normal_y;
            boundary_normal_z = old.boundary_normal_z;
            boundary_u = old.boundary_u;
            boundary_v = old.boundary_v;
            boundary_w = old.boundary_w;

            inlets = new List<Inlet>(old.inlets);

            // Re-initialize any additional members
            velocityCalculator = new VelocityCalculator();
            cellMapper = new PhysicalToCellMapper(Nx, Ny, Nz, length_x, length_y, length_z);
        }

        /// <summary>
        /// Example of how to use the PhysicalToCellMapper
        /// </summary>
        /// <param name="position">Physical position to map</param>
        /// <returns>Corresponding cell indices</returns>
        public (int i, int j, int k) MapPositionToCell(Vector3 position)
        {
            return cellMapper.GetCellIndices(position);
        }

        /// <summary>
        /// Adds a louvered inlet based on physical start/end, two tangential velocities (vT1, vT2), and total magnitude M.
        /// </summary>
        /// <param name="start">Physical starting position of the inlet.</param>
        /// <param name="end">Physical ending position of the inlet.</param>
        /// <param name="vT1">Tangential velocity component 1.</param>
        /// <param name="vT2">Tangential velocity component 2.</param>
        /// <param name="totalMagnitude">Desired total velocity magnitude M.</param>
        public void AddLouveredInlet(
            Vector3 start,
            Vector3 end,
            float vT1,
            float vT2,
            float totalMagnitude)
        {
            // Determine indices (startIndex1, startIndex2, size1, size2, wall)
            var inletParams = cellMapper.GetInletParameters(start, end);

            int startIndex1 = inletParams.startIndex1;
            int startIndex2 = inletParams.startIndex2;
            int size1 = inletParams.size1;
            int size2 = inletParams.size2;
            Wall wall = inletParams.wall;

            Console.WriteLine("Simulation completed.");

            // Print the parameters from the cell mapper
            Console.WriteLine("Inlet Parameters from Cell Mapper:");
            Console.WriteLine($"Start Index 1: {startIndex1}");
            Console.WriteLine($"Start Index 2: {startIndex2}");
            Console.WriteLine($"Size 1: {size1}");
            Console.WriteLine($"Size 2: {size2}");
            Console.WriteLine($"Wall: {wall}");

            // Get the wall's outward normal
            Vector3 normalDirection = GetNormalDirection(wall);

            // Compute velocity vector using the new method
            VelocityComponents louveredComponents = velocityCalculator.CalculateVelocityComponents(
                normalDirection,
                vT1,
                vT2,
                totalMagnitude
            );
            //Vector3 velocity = louveredComponents.TotalVelocity;
            Vector3 velocity = normalDirection * totalMagnitude;

            // Add the inlet with the computed velocity
            AddInlet(
                wall: (Inlet.Wall)wall,
                type: Inlet.InletType.Louvered,
                startIndex1: startIndex1,
                startIndex2: startIndex2,
                size1: size1,
                size2: size2,
                velocity: louveredComponents.TotalVelocity
            );
        }

        /// <summary>
        /// Determines the wall's outward normal direction as a unit vector.
        /// </summary>
        private Vector3 GetNormalDirection(Wall wall)
        {
            return wall switch
            {
                Wall.Left => Vector3.UnitX,     // +X
                Wall.Right => -Vector3.UnitX,    // -X
                Wall.Front => Vector3.UnitY,     // +Y
                Wall.Back => -Vector3.UnitY,    // -Y
                Wall.Top => -Vector3.UnitZ,    // +Z
                Wall.Bottom => Vector3.UnitZ,     // -Z
                _ => throw new ArgumentException("Invalid wall specified for normal direction."),
            };
        }

        /// <summary>
        /// Adds a standard inlet to the domain using physical positions and a total velocity magnitude.
        /// The inlet is positioned based on the physical start/end positions, and its velocity
        /// is computed as the wall’s outward normal scaled by the provided total magnitude.
        /// </summary>
        /// <param name="start">Physical starting position of the inlet.</param>
        /// <param name="end">Physical ending position of the inlet.</param>
        /// <param name="totalVelocityMagnitude">Total velocity magnitude for the inlet.</param>
        public void AddStandardInlet(Vector3 start, Vector3 end, float totalVelocityMagnitude)
        {
            // Determine indices and wall based on physical positions
            var inletParams = cellMapper.GetInletParameters(start, end);
            int startIndex1 = inletParams.startIndex1;
            int startIndex2 = inletParams.startIndex2;
            int size1 = inletParams.size1;
            int size2 = inletParams.size2;
            Wall wall = inletParams.wall;

            Console.WriteLine("Standard Inlet Parameters from Cell Mapper:");
            Console.WriteLine($"Start Index 1: {startIndex1}");
            Console.WriteLine($"Start Index 2: {startIndex2}");
            Console.WriteLine($"Size 1: {size1}");
            Console.WriteLine($"Size 2: {size2}");
            Console.WriteLine($"Wall: {wall}");

            // Get the wall's outward normal as a unit vector.
            Vector3 normalDirection = GetNormalDirection(wall);

            // Calculate the velocity vector: direction based on the wall normal scaled by total magnitude.
            Vector3 velocity = normalDirection * totalVelocityMagnitude;

            // Add the inlet using the calculated indices and velocity.
            AddInlet(
                wall: (Inlet.Wall)wall,
                type: Inlet.InletType.Standard,
                startIndex1: startIndex1,
                startIndex2: startIndex2,
                size1: size1,
                size2: size2,
                velocity: velocity
            );
        }

        // You can add more specialized inlet methods here if needed

        // Rest of the RoomVentilationDomain class remains unchanged

        /// <summary>
        /// Adds a 2D inlet to the domain, fully attached to a specified wall.
        /// </summary>
        /// <param name="wall">Wall where the inlet is located</param>
        /// <param name="type">Type of the inlet (Standard or Louvered)</param>
        /// <param name="startIndex1">Starting index in the first direction on the wall</param>
        /// <param name="startIndex2">Starting index in the second direction on the wall</param>
        /// <param name="size1">Size in the first direction (number of cells)</param>
        /// <param name="size2">Size in the second direction (number of cells)</param>
        /// <param name="velocity">Velocity vector enforced by the inlet</param>
        public void AddInlet(
            Inlet.Wall wall,
            Inlet.InletType type,
            int startIndex1,
            int startIndex2,
            int size1,
            int size2,
            Vector3 velocity)
        {
            // Validate indices based on wall
            switch (wall)
            {
                case Inlet.Wall.Left:
                case Inlet.Wall.Right:
                    if (startIndex1 < 1 || startIndex1 + size1 > Ny - 1 ||
                        startIndex2 < 1 || startIndex2 + size2 > Nz - 1)
                        throw new ArgumentException("Invalid inlet start index or size for Left/Right wall.");
                    break;
                case Inlet.Wall.Front:
                case Inlet.Wall.Back:
                    if (startIndex1 < 1 || startIndex1 + size1 > Nx - 1 ||
                        startIndex2 < 1 || startIndex2 + size2 > Nz - 1)
                        throw new ArgumentException("Invalid inlet start index or size for Front/Back wall.");
                    break;
                case Inlet.Wall.Top:
                case Inlet.Wall.Bottom:
                    if (startIndex1 < 1 || startIndex1 + size1 > Nx - 1 ||
                        startIndex2 < 1 || startIndex2 + size2 > Ny - 1)
                        throw new ArgumentException("Invalid inlet start index or size for Top/Bottom wall.");
                    break;
                default:
                    throw new ArgumentException("Invalid wall specified for inlet.");
            }

            // Add the inlet to the list
            inlets.Add(new Inlet(wall, type, startIndex1, startIndex2, size1, size2, velocity));

            // Apply boundary conditions immediately after adding the inlet
            ApplyInletBoundary(inlets[^1]); // Using C# 8.0 index from end operator

            // Log inlet addition details for verification
            Console.WriteLine($"AddInlet: Added {type} inlet on {wall} wall with size1={size1}, size2={size2}, velocity=({velocity.X}, {velocity.Y}, {velocity.Z})");
        }

        /// <summary>
        /// Applies the inlet boundary conditions based on the inlet's properties.
        /// </summary>
        /// <param name="inlet">The inlet to apply</param>
        private void ApplyInletBoundary(Inlet inlet)
        {
            switch (inlet.InletWall)
            {
                case Inlet.Wall.Left:
                    ApplyInletToWall(
                        inlet,
                        fixedIndexFunc: _ => 1, // Fixed i for Left wall
                        index1Func: a => inlet.StartIndex1 + a,
                        index2Func: b => inlet.StartIndex2 + b,
                        applyAction: (i, j, k) =>
                        {
                            boundary_u[i, j, k] = boundary_u[i, j, k] + inlet.Velocity.X;
                            boundary_v[i, j, k] = boundary_v[i, j, k] + inlet.Velocity.Y;
                            boundary_w[i, j, k] = boundary_w[i, j, k] + inlet.Velocity.Z;
                            boundary_cells[i, j, k] = 2;
                        }
                    );
                    break;

                case Inlet.Wall.Right:
                    ApplyInletToWall(
                        inlet,
                        fixedIndexFunc: _ => Nx - 2, // Fixed i for Right wall
                        index1Func: a => inlet.StartIndex1 + a,
                        index2Func: b => inlet.StartIndex2 + b,
                        applyAction: (i, j, k) =>
                        {
                            boundary_u[i, j, k] = boundary_u[i, j, k] + inlet.Velocity.X;
                            boundary_v[i, j, k] = boundary_v[i, j, k] + inlet.Velocity.Y;
                            boundary_w[i, j, k] = boundary_w[i, j, k] + inlet.Velocity.Z;
                            boundary_cells[i, j, k] = 2;
                        }
                    );
                    break;

                case Inlet.Wall.Front:
                    ApplyInletToWall(
                        inlet,
                        fixedIndexFunc: _ => 1, // Fixed j for Front wall (Y=0)
                        index1Func: a => inlet.StartIndex1 + a,
                        index2Func: b => inlet.StartIndex2 + b,
                        applyAction: (i, j, k) =>
                        {
                            boundary_u[i, j, k] = boundary_u[i, j, k] + inlet.Velocity.X;
                            boundary_v[i, j, k] = boundary_v[i, j, k] + inlet.Velocity.Y;
                            boundary_w[i, j, k] = boundary_w[i, j, k] + inlet.Velocity.Z;
                            boundary_cells[i, j, k] = 2;
                        }
                    );
                    break;

                case Inlet.Wall.Back:
                    ApplyInletToWall(
                        inlet,
                        fixedIndexFunc: _ => Ny - 2, // Fixed j for Back wall (Y=Ny-1)
                        index1Func: a => inlet.StartIndex1 + a,
                        index2Func: b => inlet.StartIndex2 + b,
                        applyAction: (i, j, k) =>
                        {
                            Console.WriteLine("Index i ={0}, j={1}, k={2}", i, j, k);
                            boundary_u[i, j, k] = boundary_u[i, j, k] + inlet.Velocity.X;
                            boundary_v[i, j, k] = boundary_v[i, j, k] + inlet.Velocity.Y;
                            boundary_w[i, j, k] = boundary_w[i, j, k] + inlet.Velocity.Z;
                            boundary_cells[i, j, k] = 2;
                        }
                    );
                    break;

                case Inlet.Wall.Top:
                    ApplyInletToWall(
                        inlet,
                        fixedIndexFunc: _ => Nz - 2, // k=Nz-2 (since indexing starts at 0)
                        index1Func: a => inlet.StartIndex1 + a, // y index
                        index2Func: b => inlet.StartIndex2 + b, // z index
                        applyAction: (i, j, k) =>
                        {
                            // Validate indices
                            bool isBoundaryUValid = i >= 0 && i < boundary_u.GetLength(0) &&
                                                    j >= 0 && j < boundary_u.GetLength(1) &&
                                                    k >= 0 && k < boundary_u.GetLength(2);
                            bool isBoundaryWValid = i >= 0 && i < boundary_w.GetLength(0) &&
                                                    j >= 0 && j < boundary_w.GetLength(1) &&
                                                    k >= 0 && k < boundary_w.GetLength(2);
                            bool isBoundaryCellsValid = i >= 0 && i < boundary_cells.GetLength(0) &&
                                                        j >= 0 && j < boundary_cells.GetLength(1) &&
                                                        k >= 0 && k < boundary_cells.GetLength(2);

                            // Log out-of-range indices
                            if (!isBoundaryUValid)
                            {
                                Console.WriteLine($"boundary_u index out of range: i={i}, j={j}, k={k}");
                            }
                            if (!isBoundaryWValid)
                            {
                                Console.WriteLine($"boundary_w index out of range: i={i}, j={j}, k={k}");
                            }
                            if (!isBoundaryCellsValid)
                            {
                                Console.WriteLine($"boundary_cells index out of range: i={i}, j={j}, k={k}");
                            }

                            // Only proceed if all indices are valid
                            if (isBoundaryUValid && isBoundaryWValid && isBoundaryCellsValid)
                            {
                                try
                                {
                                    boundary_u[i, j, k] = boundary_u[i, j, k] + inlet.Velocity.X;
                                    boundary_v[i, j, k] = boundary_v[i, j, k] + inlet.Velocity.Y;
                                    boundary_w[i, j, k] = boundary_w[i, j, k] + inlet.Velocity.Z;
                                    boundary_cells[i, j, k] = 2;
                                    //outflow_boundary_x[i, j, k] = 1;

                                    // Log the assignment
                                  //  Console.WriteLine($"Assigned inlet cell at (i={i}, j={j}, k={k}): u={boundary_u[i, j, k]},  v={boundary_v[i, j, k]}, w={boundary_w[i, j, k]}");
                                }
                                catch (IndexOutOfRangeException ex)
                                {
                                    Console.WriteLine($"Exception during assignment: {ex.Message}");
                                    throw;
                                }
                            }
                            else
                            {
                                Console.WriteLine($"Skip assignment for invalid indices: i={i}, j={j}, k={k}");
                            }
                        }
                    );
                    break;


                case Inlet.Wall.Bottom:
                    ApplyInletToWall(
                        inlet,
                        fixedIndexFunc: _ => 1, // Fixed k for Bottom wall (Z=0)
                        index1Func: a => inlet.StartIndex1 + a,
                        index2Func: b => inlet.StartIndex2 + b,
                        applyAction: (i, j, k) =>
                        {
                            // Validate indices
                            bool isBoundaryUValid = i >= 0 && i < boundary_u.GetLength(0) &&
                                                    j >= 0 && j < boundary_u.GetLength(1) &&
                                                    k >= 0 && k < boundary_u.GetLength(2);
                            bool isBoundaryWValid = i >= 0 && i < boundary_w.GetLength(0) &&
                                                    j >= 0 && j < boundary_w.GetLength(1) &&
                                                    k >= 0 && k < boundary_w.GetLength(2);
                            bool isBoundaryCellsValid = i >= 0 && i < boundary_cells.GetLength(0) &&
                                                        j >= 0 && j < boundary_cells.GetLength(1) &&
                                                        k >= 0 && k < boundary_cells.GetLength(2);

                            // Log out-of-range indices
                            if (!isBoundaryUValid)
                            {
                                Console.WriteLine($"boundary_u index out of range: i={i}, j={j}, k={k}");
                            }
                            if (!isBoundaryWValid)
                            {
                                Console.WriteLine($"boundary_w index out of range: i={i}, j={j}, k={k}");
                            }
                            if (!isBoundaryCellsValid)
                            {
                                Console.WriteLine($"boundary_cells index out of range: i={i}, j={j}, k={k}");
                            }

                            // Only proceed if all indices are valid
                            if (isBoundaryUValid && isBoundaryWValid && isBoundaryCellsValid)
                            {
                                try
                                {
                                    boundary_u[i, j, k] = boundary_u[i, j, k] + inlet.Velocity.X;
                                    boundary_v[i, j, k] = boundary_v[i, j, k] + inlet.Velocity.Y;
                                    boundary_w[i, j, k] = boundary_w[i, j, k] + inlet.Velocity.Z;
                                    boundary_cells[i, j, k] = 2;
                                }
                                catch (IndexOutOfRangeException ex)
                                {
                                    Console.WriteLine($"Exception during assignment: {ex.Message}");
                                    throw; // Rethrow after logging
                                }
                            }
                            else
                            {
                                Console.WriteLine($"Skip assignment for invalid indices: i={i}, j={j}, k={k}");
                            }
                        }
                    );
                    break;

                default:
                    throw new ArgumentException("Invalid wall specified for inlet.");
            }

            // **Count and log the number of inlet cells marked**
            int inletCellCount = 0;
            for (int i = 0; i < Nx; i++)
            {
                for (int j = 0; j < Ny; j++)
                {
                    for (int k = 0; k < Nz; k++)
                    {
                        if (boundary_cells[i, j, k] == 2)
                        {
                            inletCellCount++;
                        }
                    }
                }
            }
            Console.WriteLine($"ApplyInletBoundary: {inletCellCount} inlet cells marked.");
        }

        /// <summary>
        /// Helper method to apply inlet to a specific wall.
        /// </summary>
        /// <param name="inlet">The inlet to apply</param>
        /// <param name="fixedIndexFunc">Function to determine the fixed index based on wall</param>
        /// <param name="index1Func">Function to determine the first variable index</param>
        /// <param name="index2Func">Function to determine the second variable index</param>
        /// <param name="applyAction">Action to apply velocity and mark the cell</param>
        private void ApplyInletToWall(
            Inlet inlet,
            Func<int, int> fixedIndexFunc,
            Func<int, int> index1Func,
            Func<int, int> index2Func,
            Action<int, int, int> applyAction)
        {
            for (int a = 0; a < inlet.Size1; a++)
            {
                for (int b = 0; b < inlet.Size2; b++)
                {
                    int index1 = index1Func(a);
                    int index2 = index2Func(b);
                    int fixedIndex = fixedIndexFunc(0); // The fixed index is either i, j, or k depending on the wall

                    // Determine (i, j, k) based on the wall
                    switch (inlet.InletWall)
                    {
                        case Inlet.Wall.Left:
                        case Inlet.Wall.Right:
                            applyAction(fixedIndex, index1, index2);
                            break;
                        case Inlet.Wall.Front:
                        case Inlet.Wall.Back:
                            applyAction(index1, fixedIndex, index2);
                            break;
                        case Inlet.Wall.Top:
                        case Inlet.Wall.Bottom:
                            applyAction(index1, index2, fixedIndex);
                            break;
                        default:
                            throw new ArgumentException("Invalid wall specified for inlet.");
                    }
                }
            }
        }

        /// <summary>
        /// Reports all boundary conditions that are different from zero, including velocities.
        /// </summary>
        public void ReportNonZeroBoundaryConditions()
        {
            Console.WriteLine("=== Non-Zero Boundary Conditions Report ===");

            // Report non-zero boundary_u
            for (int i = 0; i < boundary_u.GetLength(0); i++)
            {
                for (int j = 0; j < boundary_u.GetLength(1); j++)
                {
                    for (int k = 0; k < boundary_u.GetLength(2); k++)
                    {
                        if (boundary_u[i, j, k] != 0.0)
                        {
                            string wall = DetermineWallFromBoundaryCondition("u", i, j, k);
                            double velocityMagnitude = Math.Sqrt(
                                Math.Pow(boundary_u[i, j, k], 2) +
                                Math.Pow(boundary_v[i, j, k], 2) +
                                Math.Pow(boundary_w[i, j, k], 2)
                            );
                            Console.WriteLine($"Boundary_u - Wall: {wall}, Cell: (i={i}, j={j}, k={k}), u={boundary_u[i, j, k]}, Velocity Magnitude={velocityMagnitude}");
                        }
                    }
                }
            }

            // Report non-zero boundary_v
            for (int i = 0; i < boundary_v.GetLength(0); i++)
            {
                for (int j = 0; j < boundary_v.GetLength(1); j++)
                {
                    for (int k = 0; k < boundary_v.GetLength(2); k++)
                    {
                        if (boundary_v[i, j, k] != 0.0)
                        {
                            string wall = DetermineWallFromBoundaryCondition("v", i, j, k);
                            double velocityMagnitude = Math.Sqrt(
                                Math.Pow(boundary_u[i, j, k], 2) +
                                Math.Pow(boundary_v[i, j, k], 2) +
                                Math.Pow(boundary_w[i, j, k], 2)
                            );
                            Console.WriteLine($"Boundary_v - Wall: {wall}, Cell: (i={i}, j={j}, k={k}), v={boundary_v[i, j, k]}, Velocity Magnitude={velocityMagnitude}");
                        }
                    }
                }
            }

            // Report non-zero boundary_w
            for (int i = 0; i < boundary_w.GetLength(0); i++)
            {
                for (int j = 0; j < boundary_w.GetLength(1); j++)
                {
                    for (int k = 0; k < boundary_w.GetLength(2); k++)
                    {
                        if (boundary_w[i, j, k] != 0.0)
                        {
                            string wall = DetermineWallFromBoundaryCondition("w", i, j, k);
                            double velocityMagnitude = Math.Sqrt(
                                Math.Pow(boundary_u[i, j, k], 2) +
                                Math.Pow(boundary_v[i, j, k], 2) +
                                Math.Pow(boundary_w[i, j, k], 2)
                            );
                            Console.WriteLine($"Boundary_w - Wall: {wall}, Cell: (i={i}, j={j}, k={k}), w={boundary_w[i, j, k]}, Velocity Magnitude={velocityMagnitude}");
                        }
                    }
                }
            }

            Console.WriteLine("=== End of Report ===");
        }



        /// <summary>
        /// Determines the wall based on boundary condition type and cell indices.
        /// </summary>
        /// <param name="type">Type of boundary condition ('u', 'v', 'w', or 'cell')</param>
        /// <param name="i">Cell index in x-direction</param>
        /// <param name="j">Cell index in y-direction</param>
        /// <param name="k">Cell index in z-direction</param>
        /// <returns>Wall name as string</returns>
        private string DetermineWallFromBoundaryCondition(string type, int i, int j, int k)
        {
            switch (type)
            {
                case "u":
                    if (i == 1)
                        return "Left";
                    if (i == Nx - 1)
                        return "Right";
                    break;
                case "v":
                    if (j == 1)
                        return "Front";
                    if (j == Ny - 1)
                        return "Back";
                    break;
                case "w":
                    if (k == 1)
                        return "Bottom";
                    if (k == Nz - 1)
                        return "Top";
                    break;
                case "cell":
                    // Determine wall based on fixed index in cell_cells
                    // Assuming that boundary_cells are marked as inlet cells
                    // You might need additional information to accurately determine the wall
                    // For simplicity, using the position of the cell
                    if (i == 1)
                        return "Left";
                    if (i == Nx - 1)
                        return "Right";
                    if (j == 1)
                        return "Front";
                    if (j == Ny - 1)
                        return "Back";
                    if (k == 1)
                        return "Bottom";
                    if (k == Nz - 1)
                        return "Top";
                    break;
            }
            return "Unknown";
        }
    }
}

