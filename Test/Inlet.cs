using System.Numerics;

public class Inlet
{
    public enum Wall
    {
        Left,
        Right,
        Front,
        Back,
        Top,
        Bottom
    }

public enum InletType
{
    Standard,
    Louvered,
    Radial,
    Spiral
}

    public Wall InletWall { get; set; }
    public InletType Type { get; set; }
    public int StartIndex1 { get; set; }
    public int StartIndex2 { get; set; }
    public int Size1 { get; set; }
    public int Size2 { get; set; }
    public Vector3 Velocity { get; set; }

    public Inlet(Wall wall, InletType type, int start1, int start2, int size1, int size2, Vector3 velocity)
    {
        InletWall = wall;
        Type = type;
        StartIndex1 = start1;
        StartIndex2 = start2;
        Size1 = size1;
        Size2 = size2;
        Velocity = velocity;
    }
}
