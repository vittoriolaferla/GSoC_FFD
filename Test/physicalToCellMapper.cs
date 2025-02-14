using System;
using System.Numerics;

namespace FastFluidSolver
{
    /// <summary>
    /// Enum representing the different walls of the simulation domain.
    /// </summary>
    public enum Wall
    {
        Left,
        Right,
        Front,
        Back,
        Top,
        Bottom
    }

    /// <summary>
    /// Class responsible for mapping physical positions to mesh cell indices.
    /// </summary>
    public class PhysicalToCellMapper
    {
        // Mesh dimensions (number of cells in each direction)
        public int Nx { get; }
        public int Ny { get; }
        public int Nz { get; }

        // Physical lengths of the domain (in meters)
        public double LengthX { get; }
        public double LengthY { get; }
        public double LengthZ { get; }

        // Cell sizes (in meters)
        public double Hx { get; }
        public double Hy { get; }
        public double Hz { get; }

        // Tolerance for floating-point comparisons
        private const double Tolerance = 1e-3;

        /// <summary>
        /// Initializes a new instance of the <see cref="PhysicalToCellMapper"/> class.
        /// </summary>
        /// <param name="Nx">Number of cells in the x-direction (including ghost cells).</param>
        /// <param name="Ny">Number of cells in the y-direction (including ghost cells).</param>
        /// <param name="Nz">Number of cells in the z-direction (including ghost cells).</param>
        /// <param name="length_x">Physical length of the domain in the x-direction (meters).</param>
        /// <param name="length_y">Physical length of the domain in the y-direction (meters).</param>
        /// <param name="length_z">Physical length of the domain in the z-direction (meters).</param>
        public PhysicalToCellMapper(int Nx, int Ny, int Nz, double length_x, double length_y, double length_z)
        {
            if (Nx < 3 || Ny < 3 || Nz < 3)
                throw new ArgumentException("Mesh dimensions must be at least 3 in each direction to accommodate ghost cells.");

            this.Nx = Nx;
            this.Ny = Ny;
            this.Nz = Nz;

            this.LengthX = length_x;
            this.LengthY = length_y;
            this.LengthZ = length_z;

            // Calculate cell sizes based on physical lengths and mesh counts (excluding ghost cells)
            this.Hx = length_x / (Nx - 2);
            this.Hy = length_y / (Ny - 2);
            this.Hz = length_z / (Nz - 2);
        }

        /// <summary>
        /// Transforms a physical position to corresponding cell indices in the mesh.
        /// </summary>
        /// <param name="position">Physical position as a Vector3 (meters).</param>
        /// <returns>Tuple containing cell indices (i, j, k).</returns>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when the position is outside the physical domain.</exception>
        public (int i, int j, int k) GetCellIndices(Vector3 position)
        {
            // Extract physical coordinates
            double x = position.X;
            double y = position.Y;
            double z = position.Z;

            // Validate position within the physical domain
            if (x < 0 || x > LengthX || y < 0 || y > LengthY || z < 0 || z > LengthZ)
                throw new ArgumentOutOfRangeException("Position is outside the physical domain.");

            // Calculate cell indices (1-based due to ghost cells)
            int i = (int)Math.Floor(x / Hx) + 1;
            int j = (int)Math.Floor(y / Hy) + 1;
            int k = (int)Math.Floor(z / Hz) + 1;

            // Clamp indices to ensure they fall within the valid range
            i = Math.Clamp(i, 1, Nx - 2);
            j = Math.Clamp(j, 1, Ny - 2);
            k = Math.Clamp(k, 1, Nz - 2);

            return (i, j, k);
        }

        /// <summary>
        /// Transforms cell indices back to the center physical position of the cell.
        /// </summary>
        /// <param name="i">Cell index in the x-direction.</param>
        /// <param name="j">Cell index in the y-direction.</param>
        /// <param name="k">Cell index in the z-direction.</param>
        /// <returns>Physical position as a Vector3 (meters).</returns>
        /// <exception cref="ArgumentOutOfRangeException">Thrown when the indices are out of bounds.</exception>
        public Vector3 GetPhysicalPosition(int i, int j, int k)
        {
            if (i < 1 || i >= Nx - 1 || j < 1 || j >= Ny - 1 || k < 1 || k >= Nz - 1)
                throw new ArgumentOutOfRangeException("Cell indices are out of the valid range.");

            // Calculate physical coordinates (cell centers)
            double x = (i - 0.5) * Hx;
            double y = (j - 0.5) * Hy;
            double z = (k - 0.5) * Hz;

            return new Vector3((float)x, (float)y, (float)z);
        }

        /// <summary>
        /// Determines the inlet parameters based on physical start and end positions.
        /// </summary>
        /// <param name="start">Physical starting position of the inlet.</param>
        /// <param name="end">Physical ending position of the inlet.</param>
        /// <returns>
        /// A tuple containing:
        /// - startIndex1: Starting index in the first direction on the wall.
        /// - startIndex2: Starting index in the second direction on the wall.
        /// - size1: Size in the first direction (number of cells).
        /// - size2: Size in the second direction (number of cells).
        /// - wall: The wall on which the inlet is located.
        /// </returns>
        /// <exception cref="ArgumentException">Thrown when points are not on the same wall or do not lie on a valid wall.</exception>
        public (int startIndex1, int startIndex2, int size1, int size2, Wall wall) GetInletParameters(Vector3 start, Vector3 end)
        {
            // Determine the wall for both start and end points
            Wall? startWall = DetermineWall(start);
            Wall? endWall = DetermineWall(end);

            // Validate that both points are on the same wall
            if (startWall == null || endWall == null)
                throw new ArgumentException("One or both positions do not lie on a recognized wall.");

            if (startWall != endWall)
                throw new ArgumentException("Start and end positions are not on the same wall.");

            Wall wall = startWall.Value;

            // Map physical positions to cell indices
            var startIndices = GetCellIndices(start);
            var endIndices = GetCellIndices(end);

            // Depending on the wall, determine startIndex1, startIndex2, size1, size2
            int startIndex1, startIndex2, size1, size2;

            switch (wall)
            {
                case Wall.Left:
                case Wall.Right:
                    // Variation along Y and Z axes
                    startIndex1 = Math.Min(startIndices.j, endIndices.j);
                    startIndex2 = Math.Min(startIndices.k, endIndices.k);
                    size1 = Math.Abs(endIndices.j - startIndices.j) + 1;
                    size2 = Math.Abs(endIndices.k - startIndices.k) + 1;
                    break;

                case Wall.Front:
                case Wall.Back:
                    // Variation along X and Z axes
                    startIndex1 = Math.Min(startIndices.i, endIndices.i);
                    startIndex2 = Math.Min(startIndices.k, endIndices.k);
                    size1 = Math.Abs(endIndices.i - startIndices.i) + 1;
                    size2 = Math.Abs(endIndices.k - startIndices.k) + 1;
                    break;

                case Wall.Top:
                case Wall.Bottom:
                    // Variation along X and Y axes
                    startIndex1 = Math.Min(startIndices.i, endIndices.i);
                    startIndex2 = Math.Min(startIndices.j, endIndices.j);
                    size1 = Math.Abs(endIndices.i - startIndices.i) + 1;
                    size2 = Math.Abs(endIndices.j - startIndices.j) + 1;
                    break;

                default:
                    throw new ArgumentException("Unsupported wall type.");
            }

            return (startIndex1, startIndex2, size1, size2, wall);
        }

        /// <summary>
        /// Determines the wall on which a given physical position lies.
        /// </summary>
        /// <param name="position">Physical position as a Vector3 (meters).</param>
        /// <returns>The wall if the position lies on one; otherwise, null.</returns>
        private Wall? DetermineWall(Vector3 position)
        {
            // Check for Left and Right walls (X = 0 or X = LengthX)
            if (Math.Abs(position.X) <= Tolerance)
                return Wall.Left;
            if (Math.Abs(position.X - LengthX) <= Tolerance)
                return Wall.Right;

            // Check for Front and Back walls (Y = 0 or Y = LengthY)
            if (Math.Abs(position.Y) <= Tolerance)
                return Wall.Front;
            if (Math.Abs(position.Y - LengthY) <= Tolerance)
                return Wall.Back;

            // Check for Top and Bottom walls (Z = 0 or Z = LengthZ)
            if (Math.Abs(position.Z) <= Tolerance)
                return Wall.Bottom;
            if (Math.Abs(position.Z - LengthZ) <= Tolerance)
                return Wall.Top;

            // Position does not lie on any recognized wall
            return null;
        }
    }
}
