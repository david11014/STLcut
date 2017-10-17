using System;
using System.IO;
using System.Windows.Forms;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace STLcut
{
    class Model
    {
        public List<Triangle> TrArray;
        public string Filename;
        public int TriNum;
        public BoundingBox BBox;
        public Model() {
            
            Filename = "";
            TriNum = 0;
        }
        public Model(string fn)
        {
            
            Filename = fn;
            TriNum = 0;
            LoadSTLBinnary(fn);
            BBox = new BoundingBox(ref TrArray);
        }

        public void LoadSTLBinnary()
        {
            if (Filename == "")
                return;
            else
                LoadSTLBinnary(Filename);

            return;
        }
        public void LoadSTLBinnary(string fn)
        {
            
            var stlFile = new FileStream(fn, FileMode.Open);

            if (!stlFile.CanRead)
            {
                MessageBox.Show("檔名或路徑錯誤! 無法開啟檔案!");
            }
            
            var len = (int)stlFile.Length;
            var bits = new byte[len];
            stlFile.Read(bits, 0, len);

            if (BitConverter.ToString(bits,0,5) == "solid")
            {
                MessageBox.Show("is ASCII!!");
                return;
            }

            TriNum = B2I32(bits[80], bits[81], bits[82], bits[83]);
           
            TrArray = new List<Triangle>(TriNum + 1);
            const int triByte = 50;
            for (int i = 0; i < TriNum; i++)
            {
                int head = 84 + triByte * i;
                double[] val = new double[12];

                for(int j = 0; j < 12; j++)
                {
                    val[j] = B2R64(bits[head + j * 4], bits[head + j * 4 + 1], bits[head + j * 4 + 2], bits[head + j * 4 + 3]);
                }

                Vector3D n = new Vector3D(val[0], val[1], val[2]);
                Point3D p1 = new Point3D(val[3], val[4], val[5]);
                Point3D p2 = new Point3D(val[6], val[7], val[8]);
                Point3D p3 = new Point3D(val[9], val[10], val[11]);
                Triangle T = new Triangle(p1, p2, p3, n);
                TrArray.Add(T);

            }

            if(TriNum - TrArray.Count != 0 )
            {
                MessageBox.Show("讀檔錯誤: 網格數不合");
                return;
            }           

            return;
        }
        public override string ToString()
        {
            string s = "";

            foreach(Triangle T in TrArray)
            {
                s += T.ToString() + "\n";
            }

            return s;
        }
        
        public List<Triangle> CrossTri(double face,int axis) //axis x=0 y =1 z =2
        {
            List<Triangle> TriCross = new List<Triangle>();

            if(axis == 0) //X
            {
                foreach (Triangle T in TrArray)
                {
                    if(T.MinP.x < face && T.MaxP.x > face )
                    {
                        TriCross.Add(T);
                    }
                }
            }
            else if (axis == 1) //Y
            {
                foreach (Triangle T in TrArray)
                {
                    if (T.MinP.y < face && T.MaxP.y > face)
                    {
                        TriCross.Add(T);
                    }
                }
            }
            else //Z
            {
                foreach (Triangle T in TrArray)
                {
                    if (T.MinP.z < face && T.MaxP.z > face)
                    {
                        TriCross.Add(T);
                    }
                }
            }               

            return TriCross;
        }
        
        private int B2I32(byte b1, byte b2, byte b3, byte b4)
        {
            byte[] B = { b1, b2, b3, b4 };
            return BitConverter.ToInt32(B, 0);
        }
        private double B2R64(byte b1, byte b2, byte b3, byte b4)
        {
            byte[] B = { b1, b2, b3, b4 };
            //MessageBox.Show(B[0].ToString() + " " + B[1].ToString() + " " + B[2].ToString() + " " + B[3].ToString() + " ");
            //MessageBox.Show(BitConverter.ToSingle(B, 0).ToString());
            return (double)BitConverter.ToSingle(B, 0);
        }
    }
    
    class Vector3D
    {
        public double x, y, z;

        public Vector3D(){}

        public Vector3D(double a, double b, double c)
        {
            x = a;
            y = b;
            z = c;
        }

        public static Vector3D operator +(Vector3D v1, Vector3D v2)
        {
            return new Vector3D(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
        }
        public static Vector3D operator -(Vector3D v1, Vector3D v2)
        {
            return new Vector3D(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
        }
        public static Vector3D operator *(Vector3D v1, double c)
        {
            return new Vector3D(v1.x*c, v1.y * c, v1.z * c);
        }
        public static Vector3D operator /(Vector3D v1, double c)
        {
            return new Vector3D(v1.x / c, v1.y / c, v1.z / c);
        }
     

        public static double Dot(Vector3D va, Vector3D vb)
        {
            return (double)(va.x * vb.x + va.y * vb.y + va.z * vb.z);
        }
		public static Vector3D Cross(Vector3D va, Vector3D vb)
        {
            return new Vector3D( (double)(va.y * vb.z - vb.y * va.z), (double)(va.z * vb.x - vb.z * va.x), (double)(va.x * vb.y - vb.x * va.y));
        }
        public static explicit operator Vector3D(Point3D P)
        {
            return new Vector3D(P.x, P.y, P.z);
        }
        public double length()
        {
            return Math.Sqrt(x * x + y * y + z * z);
        }

        public Vector3D unit()
        {
            return new Vector3D(x/length(), y / length(), z / length());
        }

        public override string ToString()
        {
            string s = "";

            s += x.ToString() + " " + y.ToString() + " " + z.ToString();

            return s;
        }

    }

    class Point3D
    {
        public double x, y, z;

        public Point3D()
        { }

        public Point3D(double a, double b, double c)
        {
            x = a;
            y = b;
            z = c;
        }
      
        public static Point3D operator +(Point3D v1, Point3D v2)
        {
            return new Point3D(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
        }
        public static Point3D operator -(Point3D v1, Point3D v2)
        {
            return new Point3D(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
        }

        public static Point3D operator *(Point3D v1, double c)
        {
            return new Point3D(v1.x * c, v1.y * c, v1.z * c);
        }
        public static Point3D operator /(Point3D v1, double c)
        {
            return new Point3D(v1.x / c, v1.y / c, v1.z / c);
        }

        public static explicit operator Point3D(Vector3D V)
        {
            return new Point3D(V.x, V.y, V.z);
        }

        public override string ToString()
        {
            string s = "";

            s += x.ToString() + " " + y.ToString() + " " + z.ToString();

            return s;
        }

    }

    class Triangle
    {
        public Vector3D n;
        public Point3D[] p = new Point3D[3];
        public Point3D MaxP, MinP;


        //public Triangle() { }
        public Triangle(Point3D p1, Point3D p2, Point3D p3, Vector3D N)
        {
            p[0] = p1;
            p[1] = p2;
            p[2] = p3;
            n = N;

            steMaxMin();
        }

        public override string ToString()
        {
            string s = "";

            s += "n: " + n.ToString() + "\n";
            s += "p1: " + p[0].ToString() + "\n";
            s += "p2: " + p[1].ToString() + "\n";
            s += "p3: " + p[2].ToString();

            return s;
        }

        public int CompareTo(Triangle comparePart)
        {
            // A null value means that this object is greater.
            if (comparePart == null)
                return 1;

            else
                return this.p[0].y.CompareTo(comparePart.p[0].y);
        }
        private void steMaxMin()
        {
            double maxX = Math.Max(Math.Max(p[0].x, p[1].x), p[2].x);
            double maxY = Math.Max(Math.Max(p[0].y, p[1].y), p[2].y);
            double maxZ = Math.Max(Math.Max(p[0].z, p[1].z), p[2].z);

            double minX = Math.Min(Math.Min(p[0].x, p[1].x), p[2].x);
            double minY = Math.Min(Math.Min(p[0].y, p[1].y), p[2].y);
            double minZ = Math.Min(Math.Min(p[0].z, p[1].z), p[2].z);

            MaxP = new Point3D(maxX, maxY, maxZ);
            MinP = new Point3D(minX, minY, minZ);

        }

       
    }

    class BoundingBox
    {
        public Point3D MaxP, MinP;
        public BoundingBox() { }
        public BoundingBox(ref List<Triangle> Tri)
        {
            SetBoundingBox(ref Tri);
        }
        
        public void SetBoundingBox(ref List<Triangle> Tri)
        {
            double MaxX = -Double.MaxValue, MaxY = -Double.MaxValue, MaxZ = -Double.MaxValue;
            double MinX = Double.MaxValue, MinY = Double.MaxValue, MinZ = Double.MaxValue;

            foreach(Triangle T in Tri)
            {
                if (T.MaxP.x > MaxX)
                    MaxX = T.MaxP.x;

                if (T.MaxP.y > MaxY)
                    MaxY = T.MaxP.y;

                if (T.MaxP.z > MaxZ)
                    MaxZ = T.MaxP.z;

                
                if (T.MinP.x < MinX)
                    MinX = T.MinP.x;

                if (T.MinP.y < MinY)
                    MinY = T.MinP.y;

                if (T.MinP.z < MinZ)
                    MinZ = T.MinP.z;
            }

            MaxP = new Point3D(MaxX, MaxY, MaxZ);
            MinP = new Point3D(MinX, MinY, MinZ);
        }

    }


    
}
