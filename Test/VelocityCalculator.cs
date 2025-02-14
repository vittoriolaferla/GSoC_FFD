using System;
using System.Numerics;

namespace FastFluidSolver
{
    /// <summary>
    /// Stores the results of the velocity calculation (components and total).
    /// </summary>
    public struct VelocityComponents
    {
        public Vector3 TotalVelocity { get; set; }
    }

    /// <summary>
    /// Responsible for calculating the velocity vector from a normal direction,
    /// two tangential velocities, and a total desired magnitude.
    /// </summary>
    public class VelocityCalculator
    {
        /// <summary>
        /// Calculates the 3D velocity vector given:
        ///   - A normal direction (unit vector).
        ///   - Two tangential velocities (vT1, vT2) along orthonormal directions perpendicular to the normal.
        ///   - A total velocity magnitude M.
        /// </summary>
        /// <param name="normalDirection">
        ///   A unit vector normal to the wall. 
        ///   (Must be normalized before calling this method.)
        /// </param>
        /// <param name="vT1">Tangential velocity component along t1.</param>
        /// <param name="vT2">Tangential velocity component along t2.</param>
        /// <param name="totalMagnitude">Total magnitude M of the final velocity vector.</param>
        /// <returns>A <see cref="VelocityComponents"/> with the computed 3D velocity.</returns>
        /// <exception cref="ArgumentException">
        ///   Thrown if vT1^2 + vT2^2 > M^2 and we cannot compute a real normal component.
        /// </exception>
        public VelocityComponents CalculateVelocityComponents(
            Vector3 normalDirection, 
            float vT1, 
            float vT2, 
            float totalMagnitude)
        {
            // Ensure the provided normal direction is indeed normalized (or nearly so).
            float nLength = normalDirection.Length();
            if (Math.Abs(nLength - 1.0f) > 1e-4f)
            {
                // If not normalized, normalize it (or throw).
                normalDirection = Vector3.Normalize(normalDirection);
            }

            // Find two orthonormal vectors (t1, t2) perpendicular to "normalDirection".
            // A robust approach: pick any vector not parallel to normalDirection, then use cross products.
            Vector3 t1 = FindAnyPerpendicular(normalDirection);
            t1 = Vector3.Normalize(t1);
            Vector3 t2 = Vector3.Cross(normalDirection, t1);
            t2 = Vector3.Normalize(t2);

            // Compute tangential magnitude and check feasibility.
            float tangentialSq = vT1 * vT1 + vT2 * vT2;
            float totalSq = totalMagnitude * totalMagnitude;
            if (tangentialSq > totalSq)
            {
                // Depending on your model, you can clamp or throw an exception:
                throw new ArgumentException("Tangential speeds exceed total velocity magnitude. Cannot compute real normal component.");
            }

            // Normal component magnitude
            float vNormal = (float)Math.Sqrt(totalSq - tangentialSq);

            // Build final velocity
            Vector3 velocity = vNormal * normalDirection
                             + vT1 * t1
                             + vT2 * t2;

            return new VelocityComponents
            {
                TotalVelocity = velocity
            };
        }

        /// <summary>
        /// Finds any vector that is perpendicular to the provided normal direction.
        /// </summary>
        private Vector3 FindAnyPerpendicular(Vector3 n)
        {
            // If n is near the X axis, pick a Y-based vector, otherwise pick an X-based vector, etc.
            // This ensures we don't pick a parallel direction.
            if (Math.Abs(n.X) < 0.5f)
                return new Vector3(1, 0, 0) - Vector3.Dot(new Vector3(1, 0, 0), n) * n;
            else
                return new Vector3(0, 1, 0) - Vector3.Dot(new Vector3(0, 1, 0), n) * n;
        }
    }
}
