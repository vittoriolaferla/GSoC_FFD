using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

/*
 * PostProcessor.cs
 * Copyright 2016 Lukas Bystricky <lb13f@my.fsu.edu>
 *
 * This work is licensed under the GNU GPL license version 2 or later.
 */

namespace FastFluidSolver
{
    /// <summary>
    /// Class for PostProcessing data calculated by the FFD solver in FluidSolver.cs.
    /// Exports data to VTK files for further analysis. VTK files
    /// are used often for visualization of scientific data. These files can be read
    /// by Visit, Paraview and other applications. A toolkit also exists in C++ for analyzing
    /// this data. More info on VTK is avaiable here: http://www.vtk.org.
    /// </summary>
    public class PostProcessor
    {
        FluidSolver fs;
        Domain omega;
        DataExtractor de;

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="fs">FluidSolver containing solution</param>
        /// <param name="omega">Domain containing geometry information and 
        /// exact solutions (if known)</param>
        public PostProcessor(FluidSolver fs, Domain omega)
        {
            this.fs = fs;
            this.omega = omega;

            de = new DataExtractor(omega, fs);
        }

        /// <summary>
        /// Exports data to a VTK file. Data can be cell centred or vertex centred. 
        /// </summary>
        /// <param name="fname">file name</param>
        /// <param name="time">time</param>
        /// <param name="cell_centred">cell centred flag</param>
        /// <remarks>The file format guide for VTK is avaliable here:
        /// http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf </remarks>
        public void export_data_vtk(String fname, double time, bool cell_centred)
        {
            int Nx, Ny, Nz;

            double hx = omega.hx;
            double hy = omega.hy;
            double hz = omega.hz;

            double[,,] u_interp;
            double[,,] v_interp;
            double[,,] w_interp;
            double[,,] p_interp;
            double[,,] div_interp;

            double[,,] err_u;
            double[,,] err_v;
            double[,,] err_w;
            double[,,] err_p;

            interpolate_errors(out err_u, out err_v, out err_w, out err_p, time);

            if (cell_centred)
            {
                Nx = omega.Nx - 2;
                Ny = omega.Ny - 2;
                Nz = omega.Nz - 2;

                interpolate_to_cell_centre_grid(out u_interp, out v_interp, out w_interp, out p_interp, out div_interp);
            }
            else
            {

                Nx = omega.Nx - 1;
                Ny = omega.Ny - 1;
                Nz = omega.Nz - 1;

                interpolate_to_vertex_grid(out u_interp, out v_interp, out w_interp, out p_interp);
                div_interp = new double[1, 1, 1]; //have to assign div_interp to avoid compiling errors
            }

            using (StreamWriter sw = new StreamWriter(fname))
            {
                sw.WriteLine("# vtk DataFile Version 3.0");
                sw.WriteLine("Fast Fluid Dynamics data\n");
                sw.WriteLine("ASCII");
                sw.WriteLine("DATASET RECTILINEAR_GRID");
                sw.WriteLine("FIELD FieldData 1");
                sw.WriteLine("TIME 1 1 double");
                sw.WriteLine("{0}", time);
                sw.WriteLine("DIMENSIONS {0} {1} {2}", omega.Nx - 1, omega.Ny - 1, omega.Nz - 1);
                sw.WriteLine("X_COORDINATES {0} double", omega.Nx - 1);

                for (int i = 0; i < omega.Nx - 1; i++)
                {
                    sw.WriteLine("{0}", hx * i);
                }
                sw.WriteLine("Y_COORDINATES {0} double", omega.Ny - 1);
                for (int i = 0; i < omega.Ny - 1; i++)
                {
                    sw.WriteLine("{0}", hy * i);
                }
                sw.WriteLine("Z_COORDINATES {0} double", omega.Nz - 1);
                for (int i = 0; i < omega.Nz - 1; i++)
                {
                    sw.WriteLine("{0}", hz * i);
                }

                if (cell_centred)
                {
                    sw.WriteLine("CELL_DATA {0}", Nx * Ny * Nz);
                }
                else
                {
                    sw.WriteLine("POINT_DATA {0}", Nx * Ny * Nz);
                }

                sw.WriteLine("VECTORS velocity double");
                for (int k = 0; k < Nz; k++)
                {
                    for (int j = 0; j < Ny; j++)
                    {
                        for (int i = 0; i < Nx; i++)
                        {
                            sw.WriteLine("{0} {1} {2}", u_interp[i, j, k], v_interp[i, j, k], w_interp[i, j, k]);
                        }
                    }
                }

                sw.WriteLine("SCALARS pressure double {0}", 1);
                sw.WriteLine("LOOKUP_TABLE default");
                for (int k = 0; k < Nz; k++)
                {
                    for (int j = 0; j < Ny; j++)
                    {
                        for (int i = 0; i < Nx; i++)
                        {
                            sw.Write("{0} ", p_interp[i, j, k]);
                        }
                    }
                }

                if (cell_centred)
                {
                    sw.WriteLine("SCALARS div_u double {0}", 1);
                    sw.WriteLine("LOOKUP_TABLE default");
                    for (int k = 0; k < Nz; k++)
                    {
                        for (int j = 0; j < Ny; j++)
                        {
                            for (int i = 0; i < Nx; i++)
                            {
                                sw.Write("{0} ", div_interp[i, j, k]);
                            }
                        }
                    }
                }
                else
                {
                    sw.WriteLine("SCALARS err_u double {0}", 1);
                    sw.WriteLine("LOOKUP_TABLE default");
                    for (int k = 0; k < Nz; k++)
                    {
                        for (int j = 0; j < Ny; j++)
                        {
                            for (int i = 0; i < Nx; i++)
                            {
                                sw.Write("{0} ", err_u[i, j, k]);
                            }
                        }
                    }

                    sw.WriteLine("SCALARS err_v double {0}", 1);
                    sw.WriteLine("LOOKUP_TABLE default");
                    for (int k = 0; k < Nz; k++)
                    {
                        for (int j = 0; j < Ny; j++)
                        {
                            for (int i = 0; i < Nx; i++)
                            {
                                sw.Write("{0} ", err_v[i, j, k]);
                            }
                        }
                    }

                    sw.WriteLine("SCALARS err_w double {0}", 1);
                    sw.WriteLine("LOOKUP_TABLE default");
                    for (int k = 0; k < Nz; k++)
                    {
                        for (int j = 0; j < Ny; j++)
                        {
                            for (int i = 0; i < Nx; i++)
                            {
                                sw.Write("{0} ", err_w[i, j, k]);
                            }
                        }
                    }

                    sw.WriteLine("SCALARS err_p double {0}", 1);
                    sw.WriteLine("LOOKUP_TABLE default");
                    for (int k = 0; k < Nz; k++)
                    {
                        for (int j = 0; j < Ny; j++)
                        {
                            for (int i = 0; i < Nx; i++)
                            {
                                sw.Write("{0} ", err_p[i, j, k]);
                            }
                        }
                    }
                }
            }
        }


        public void export_intermediate_vtk(string fname, double time, Boolean cell_centred)
        {
            int Nx = fs.u_intermediate.GetLength(0);
            int Ny = fs.u_intermediate.GetLength(1);
            int Nz = fs.u_intermediate.GetLength(2);

            double hx = omega.hx;
            double hy = omega.hy;
            double hz = omega.hz;

            if (cell_centred)
            {
                Nx = omega.Nx - 2;
                Ny = omega.Ny - 2;
                Nz = omega.Nz - 2;

            }
            else
            {

                Nx = omega.Nx - 1;
                Ny = omega.Ny - 1;
                Nz = omega.Nz - 1;
            }

            using (StreamWriter sw = new StreamWriter(fname))
            {
                sw.WriteLine("# vtk DataFile Version 3.0");
                sw.WriteLine("FFD intermediate velocity\n");
                sw.WriteLine("ASCII");
                sw.WriteLine("DATASET RECTILINEAR_GRID");
                sw.WriteLine("FIELD FieldData 1");
                sw.WriteLine("TIME 1 1 double");
                sw.WriteLine("{0}", time);

                // In a rectilinear grid, specify the grid dimensions
                // If your solver's final velocity is (Nx,Ny,Nz), you might need Nx+1 etc.
                // Adjust as needed so it matches how you do it in export_data_vtk(...)
                sw.WriteLine("DIMENSIONS {0} {1} {2}", Nx, Ny, Nz);

                // X_COORDINATES
                sw.WriteLine("X_COORDINATES {0} double", Nx);
                for (int i = 0; i < Nx; i++) sw.WriteLine("{0}", i * hx);

                // Y_COORDINATES
                sw.WriteLine("Y_COORDINATES {0} double", Ny);
                for (int j = 0; j < Ny; j++) sw.WriteLine("{0}", j * hy);

                // Z_COORDINATES
                sw.WriteLine("Z_COORDINATES {0} double", Nz);
                for (int k = 0; k < Nz; k++) sw.WriteLine("{0}", k * hz);

                // For a rectilinear grid, you typically do POINT_DATA or CELL_DATA
                // depending on how you interpret Nx,Ny,Nz above. 
                // Let's do POINT_DATA here:
                sw.WriteLine("POINT_DATA {0}", Nx * Ny * Nz);

                // Now write the velocity vector
                sw.WriteLine("VECTORS velocity_intermediate double");
                for (int k = 0; k < Nz; k++)
                {
                    for (int j = 0; j < Ny; j++)
                    {
                        for (int i = 0; i < Nx; i++)
                        {
                            // At each point (i,j,k), write (u_intermediate, v_intermediate, w_intermediate)
                            sw.WriteLine("{0} {1} {2}",
                                fs.u_intermediate[i, j, k],
                                fs.v_intermediate[i, j, k],
                                fs.w_intermediate[i, j, k]);
                        }
                    }
                }
            }
        }



        /// <summary>
        /// Exports geometry, including obstacle cells (excluding ghost cells), boundary cells,
        /// outflow cells, and inlet cells to a VTK file.
        /// </summary>
        /// <param name="fname">File name for the VTK output.</param>
        /// <param name="time">Current simulation time.</param>
        public void export_geometry_vtk(String fname, double time)
        {
            double hx = omega.length_x / (omega.Nx - 2);
            double hy = omega.length_y / (omega.Ny - 2);
            double hz = omega.length_z / (omega.Nz - 2);

            using (StreamWriter sw = new StreamWriter(fname))
            {
                // VTK Header
                sw.WriteLine("# vtk DataFile Version 3.0");
                sw.WriteLine("Fast Fluid Dynamics geometry");
                sw.WriteLine("ASCII");
                sw.WriteLine("DATASET RECTILINEAR_GRID");

                // DIMENSIONS
                sw.WriteLine("DIMENSIONS {0} {1} {2}", omega.Nx + 1, omega.Ny + 1, omega.Nz + 1);

                // X_COORDINATES
                sw.WriteLine("X_COORDINATES {0} double", omega.Nx + 1);
                for (int i = 0; i < omega.Nx + 1; i++)
                {
                    sw.Write("{0} ", hx * (i - 1));
                }
                sw.WriteLine(); // Ensure newline after coordinates

                // Y_COORDINATES
                sw.WriteLine("Y_COORDINATES {0} double", omega.Ny + 1);
                for (int i = 0; i < omega.Ny + 1; i++)
                {
                    sw.Write("{0} ", hy * (i - 1));
                }
                sw.WriteLine(); // Ensure newline after coordinates

                // Z_COORDINATES
                sw.WriteLine("Z_COORDINATES {0} double", omega.Nz + 1);
                for (int i = 0; i < omega.Nz + 1; i++)
                {
                    sw.Write("{0} ", hz * (i - 1));
                }
                sw.WriteLine(); // Ensure newline after coordinates

                // CELL_DATA
                sw.WriteLine("CELL_DATA {0}", omega.Nx * omega.Ny * omega.Nz);

                // 1. boundary_cells
                sw.WriteLine("SCALARS boundary_cells double 1");
                sw.WriteLine("LOOKUP_TABLE default");
                for (int k = 0; k < omega.Nz; k++)
                {
                    for (int j = 0; j < omega.Ny; j++)
                    {
                        for (int i = 0; i < omega.Nx; i++)
                        {
                            sw.Write("{0} ", (double)omega.boundary_cells[i, j, k]);
                        }
                    }
                }
                sw.WriteLine(); // Newline after SCALARS

                // 2. obstacle_cells
                sw.WriteLine("SCALARS obstacle_cells double 1");
                sw.WriteLine("LOOKUP_TABLE default");
                for (int k = 0; k < omega.Nz; k++)
                {
                    for (int j = 0; j < omega.Ny; j++)
                    {
                        for (int i = 0; i < omega.Nx; i++)
                        {
                            // Only internal cells have obstacle data
                            if (i > 0 && i < omega.Nx - 1
                                && j > 0 && j < omega.Ny - 1
                                && k > 0 && k < omega.Nz - 1)
                            {
                                sw.Write("{0} ", (double)omega.obstacle_cells[i, j, k]);
                            }
                            else
                            {
                                sw.Write("0 ");
                            }
                        }
                    }
                }
                sw.WriteLine(); // Newline after SCALARS

                // 3. boundary_normal_x
                sw.WriteLine("SCALARS boundary_normal_x double 1");
                sw.WriteLine("LOOKUP_TABLE default");
                for (int k = 0; k < omega.Nz; k++)
                {
                    for (int j = 0; j < omega.Ny; j++)
                    {
                        for (int i = 0; i < omega.Nx; i++)
                        {
                            sw.Write("{0} ", (double)omega.boundary_normal_x[i, j, k]);
                        }
                    }
                }
                sw.WriteLine(); // Newline after SCALARS

                // 4. boundary_normal_y
                sw.WriteLine("SCALARS boundary_normal_y double 1");
                sw.WriteLine("LOOKUP_TABLE default");
                for (int k = 0; k < omega.Nz; k++)
                {
                    for (int j = 0; j < omega.Ny; j++)
                    {
                        for (int i = 0; i < omega.Nx; i++)
                        {
                            sw.Write("{0} ", (double)omega.boundary_normal_y[i, j, k]);
                        }
                    }
                }
                sw.WriteLine(); // Newline after SCALARS

                // 5. boundary_normal_z
                sw.WriteLine("SCALARS boundary_normal_z double 1");
                sw.WriteLine("LOOKUP_TABLE default");
                for (int k = 0; k < omega.Nz; k++)
                {
                    for (int j = 0; j < omega.Ny; j++)
                    {
                        for (int i = 0; i < omega.Nx; i++)
                        {
                            sw.Write("{0} ", (double)omega.boundary_normal_z[i, j, k]);
                        }
                    }
                }
                sw.WriteLine(); // Newline after SCALARS

                // 6. outflow_boundary_x
                sw.WriteLine("SCALARS outflow_boundary_x double 1");
                sw.WriteLine("LOOKUP_TABLE default");
                for (int k = 0; k < omega.Nz; k++)
                {
                    for (int j = 0; j < omega.Ny; j++)
                    {
                        for (int i = 0; i < omega.Nx; i++)
                        {
                            sw.Write("{0} ", (double)omega.outflow_boundary_x[i, j, k]);
                        }
                    }
                }
                sw.WriteLine(); // Newline after SCALARS

                // 7. outflow_boundary_y
                sw.WriteLine("SCALARS outflow_boundary_y double 1");
                sw.WriteLine("LOOKUP_TABLE default");
                for (int k = 0; k < omega.Nz; k++)
                {
                    for (int j = 0; j < omega.Ny; j++)
                    {
                        for (int i = 0; i < omega.Nx; i++)
                        {
                            sw.Write("{0} ", (double)omega.outflow_boundary_y[i, j, k]);
                        }
                    }
                }
                sw.WriteLine(); // Newline after SCALARS

                // 8. outflow_boundary_z
                sw.WriteLine("SCALARS outflow_boundary_z double 1");
                sw.WriteLine("LOOKUP_TABLE default");
                for (int k = 0; k < omega.Nz; k++)
                {
                    for (int j = 0; j < omega.Ny; j++)
                    {
                        for (int i = 0; i < omega.Nx; i++)
                        {
                            sw.Write("{0} ", (double)omega.outflow_boundary_z[i, j, k]);
                        }
                    }
                }
                sw.WriteLine(); // Newline after SCALARS

                // 9. inlet_cells
                sw.WriteLine("SCALARS inlet_cells double 1");
                sw.WriteLine("LOOKUP_TABLE default");
                for (int k = 0; k < omega.Nz; k++)
                {
                    for (int j = 0; j < omega.Ny; j++)
                    {
                        for (int i = 0; i < omega.Nx; i++)
                        {
                            // Mark inlet cells (boundary_cells == 2) as 1.0, else 0.0
                            double inlet_flag = (omega.boundary_cells[i, j, k] == 2) ? 1.0 : 0.0;
                            sw.Write("{0} ", inlet_flag);
                        }
                    }
                }
                sw.WriteLine(); // Newline after SCALARS
            }
        }



        /// <summary>
        /// Exports an uninterpolated pressure field to a VTK file
        /// </summary>
        /// <param name="fname">file name</param>
        /// <param name="time">time</param>
        public void export_uninterpolated_vtk(String fname, double time)
        {
            int Nx = fs.p.GetLength(0);
            int Ny = fs.p.GetLength(1);
            int Nz = fs.p.GetLength(2);

            double hx = omega.hx;
            double hy = omega.hy;
            double hz = omega.hz;

            using (StreamWriter sw = new StreamWriter(fname))
            {
                sw.WriteLine("# vtk DataFile Version 3.0");
                sw.WriteLine("Fast Fluid Dynamics data\n");
                sw.WriteLine("ASCII");
                sw.WriteLine("DATASET RECTILINEAR_GRID");
                sw.WriteLine("FIELD FieldData 1");
                sw.WriteLine("TIME 1 1 double");
                sw.WriteLine("{0}", time);
                sw.WriteLine("DIMENSIONS {0} {1} {2}", omega.Nx + 1, omega.Ny + 1, omega.Nz + 1);
                sw.WriteLine("X_COORDINATES {0} double", omega.Nx + 1);

                for (int i = 0; i < omega.Nx + 1; i++)
                {
                    sw.WriteLine("{0}", (i - 1) * hx);
                }
                sw.WriteLine("Y_COORDINATES {0} double", omega.Ny + 1);
                for (int i = 0; i < omega.Ny + 1; i++)
                {
                    sw.WriteLine("{0}", (i - 1) * hy);
                }
                sw.WriteLine("Z_COORDINATES {0} double", omega.Nz + 1);
                for (int i = 0; i < omega.Nz + 1; i++)
                {
                    sw.WriteLine("{0}", (i - 1) * hz);
                }

                sw.WriteLine("CELL_DATA {0}", Nx * Ny * Nz);
                sw.WriteLine("SCALARS p_uninterp double {0}", 1);
                sw.WriteLine("LOOKUP_TABLE default");
                for (int k = 0; k < Nz; k++)
                {
                    for (int j = 0; j < Ny; j++)
                    {
                        for (int i = 0; i < Nx; i++)
                        {
                            sw.Write("{0} ", fs.p[i, j, k]);
                        }
                    }
                }
            }
        }


        /// <summary>
        /// Interpolated velocity components, pressure and divergence field to cell centred grid.
        /// </summary>
        /// <param name="u_interp">interpolated x component of velocity</param>
        /// <param name="v_interp">interpolated y component of velocity</param>
        /// <param name="w_interp">interpolated z component of velocity</param>
        /// <param name="p_interp">interpolated pressure</param>
        /// <param name="div_interp">interpolated divergence</param>
        private void interpolate_to_cell_centre_grid(out double[,,] u_interp, out double[,,] v_interp,
                    out double[,,] w_interp, out double[,,] p_interp, out double[,,] div_interp)
        {
            int Nx = omega.Nx - 2;
            int Ny = omega.Ny - 2;
            int Nz = omega.Nz - 2;

            u_interp = new double[Nx, Ny, Nz];
            v_interp = new double[Nx, Ny, Nz];
            w_interp = new double[Nx, Ny, Nz];
            p_interp = new double[Nx, Ny, Nz];
            div_interp = new double[Nx, Ny, Nz];

            for (int i = 1; i < Nx + 1; i++)
            {
                for (int j = 1; j < Ny + 1; j++)
                {
                    for (int k = 1; k < Nz + 1; k++)
                    {
                        //if (omega.obstacle_cells[i, j, k] == 0)
                        {
                            double[] velocity_interp = de.get_velocity((i - 0.5) * omega.hx, (j - 0.5) * omega.hy, (k - 0.5) * omega.hz);
                            u_interp[i - 1, j - 1, k - 1] = velocity_interp[0];
                            v_interp[i - 1, j - 1, k - 1] = velocity_interp[1];
                            w_interp[i - 1, j - 1, k - 1] = velocity_interp[2];

                            p_interp[i - 1, j - 1, k - 1] = de.get_pressure((i - 0.5) * omega.hx, (j - 0.5) * omega.hy, (k - 0.5) * omega.hz);

                            div_interp[i - 1, j - 1, k - 1] = (fs.u[i, j, k] - fs.u[i - 1, j, k]) / omega.hx + (fs.v[i, j, k] - fs.v[i, j - 1, k]) / omega.hy +
                                    (fs.w[i, j, k] - fs.w[i, j, k - 1]) / omega.hz;
                        }
                    }
                }
            }
        }


        /// <summary>
        /// Interpolates the velocity and pressure from the staggered grid to the vertex grid.
        /// </summary>
        /// <param name="u_interp">Interpolated u-component velocity on vertex grid.</param>
        /// <param name="v_interp">Interpolated v-component velocity on vertex grid.</param>
        /// <param name="w_interp">Interpolated w-component velocity on vertex grid.</param>
        /// <param name="p_interp">Interpolated pressure on vertex grid.</param>
        private void interpolate_to_vertex_grid(out double[,,] u_interp, out double[,,] v_interp,
                    out double[,,] w_interp, out double[,,] p_interp)
        {
            int Nx = omega.Nx - 1;
            int Ny = omega.Ny - 1;
            int Nz = omega.Nz - 1;

            u_interp = new double[Nx, Ny, Nz];
            v_interp = new double[Nx, Ny, Nz];
            w_interp = new double[Nx, Ny, Nz];
            p_interp = new double[Nx, Ny, Nz];

            // Definisci la cella di interesse per il logging
            const int target_i = 20;
            const int target_j = 20;
            const int target_k = 48;

            for (int i = 0; i < Nx; i++)
            {
                for (int j = 0; j < Ny; j++)
                {
                    for (int k = 0; k < Nz; k++)
                    {
                        // Calcola le coordinate reali del vertice
                        double x = i * omega.hx;
                        double y = j * omega.hy;
                        double z = k * omega.hz;

                        // Ottieni le velocità interpolate
                        double[] velocity_interp = de.get_velocity(x, y, z);
                        u_interp[i, j, k] = velocity_interp[0];
                        v_interp[i, j, k] = velocity_interp[1];
                        w_interp[i, j, k] = velocity_interp[2];

                        // Ottieni la pressione interpolata
                        p_interp[i, j, k] = de.get_pressure(x, y, z);

                        // Log delle assegnazioni per la cella di interesse
                        if (i == target_i && j == target_j && k == target_k)
                        {
                            Console.WriteLine($"Interpolated velocity at (i={i}, j={j}, k={k}): u={u_interp[i, j, k]}, v={v_interp[i, j, k]}, w={w_interp[i, j, k]}");
                        }
                    }
                }
            }
        }



        /// <summary>
        /// Interpolates the errors in u, v, w and p on a Domain on which an exact solution is known.
        /// </summary>
        /// <param name="err_u">error in x component of velocity</param>
        /// <param name="err_v">error in ycomponent of velocity</param>
        /// <param name="err_w">error in z component of velocity</param>
        /// <param name="err_p">error in pressure</param>
        /// <param name="t">time</param>
        private void interpolate_errors(out double[,,] err_u, out double[,,] err_v,
            out double[,,] err_w, out double[,,] err_p, double t)
        {
            int Nx = omega.Nx - 1;
            int Ny = omega.Ny - 1;
            int Nz = omega.Nz - 1;


            err_u = new double[Nx, Ny, Nz];
            err_v = new double[Nx, Ny, Nz];
            err_w = new double[Nx, Ny, Nz];
            err_p = new double[Nx, Ny, Nz];

            for (int i = 0; i < Nx; i++)
            {
                for (int j = 0; j < Ny; j++)
                {
                    for (int k = 0; k < Nz; k++)
                    {
                        double x = i * omega.hx;
                        double y = j * omega.hy;
                        double z = k * omega.hz;

                        double[] coordinate = new double[] { x, y, z };

                        double[] velocity_interp = de.get_velocity(x, y, z);

                        double u_exact, v_exact, w_exact, p_exact;

                        omega.exact_solution(coordinate, fs.nu, t, out u_exact, out v_exact,
                                out w_exact, out p_exact);

                        err_u[i, j, k] = velocity_interp[0] - u_exact;

                        err_v[i, j, k] = velocity_interp[1] - v_exact;

                        err_w[i, j, k] = velocity_interp[2] - w_exact;

                        err_p[i, j, k] = de.get_pressure(x, y, z) - p_exact;
                    }
                }
            }
        }
    }
}

