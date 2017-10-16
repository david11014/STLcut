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
        public Model() {
            
            Filename = "";
            TriNum = 0;
        }
        public Model(string fn)
        {
            
            Filename = fn;
            TriNum = 0;
            LoadSTLBinnary(fn);
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
                        
            TriNum = B2I32(bits[80], bits[81], bits[82], bits[83]);
            //MessageBox.Show(len.ToString() + " " + TriNum.ToString());
            TrArray = new List<Triangle>(TriNum + 1);
            const int triByte = 50;
            for (int i = 0; i < TriNum; i++)
            {
                int head = 84 + triByte * i;
                float[] val = new float[12];

                for(int j = 0; j < 12; j++)
                {
                    val[j] = B2R32(bits[head + j * 4], bits[head + j * 4 + 1], bits[head + j * 4 + 2], bits[head + j * 4 + 3]);
                }

                Vector3D n = new Vector3D(val[0], val[1], val[2]);
                Point3D p1 = new Point3D(val[3], val[4], val[5]);
                Point3D p2 = new Point3D(val[6], val[7], val[8]);
                Point3D p3 = new Point3D(val[9], val[10], val[11]);
                Triangle T = new Triangle(p1, p2, p3, n);
                TrArray.Add(T);

            }

            return;
        }

        private int B2I32(byte b1, byte b2, byte b3, byte b4)
        {
            byte[] B = { b1, b2, b3, b4 };
            return BitConverter.ToInt32(B, 0);
        }
        private float B2R32(byte b1, byte b2, byte b3, byte b4)
        {
            byte[] B = { b1, b2, b3, b4 };
            //MessageBox.Show(B[0].ToString() + " " + B[1].ToString() + " " + B[2].ToString() + " " + B[3].ToString() + " ");
            //MessageBox.Show(BitConverter.ToSingle(B, 0).ToString());
            return BitConverter.ToSingle(B, 0);
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
        public Vector3D(float a, float b, float c)
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
            return va.x * vb.x + va.y * vb.y + va.z * vb.z;
        }
		public static Vector3D Cross(Vector3D va, Vector3D vb)
        {
            return new Vector3D( va.y * vb.z - vb.y * va.z, va.z * vb.x - vb.z * va.x, va.x * vb.y - vb.x * va.y);
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
        public Point3D(float a, float b, float c)
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

    }

    class Triangle
    {
        public Vector3D n;
        public Point3D[] p = new Point3D[3];


        public Triangle() { }
        public Triangle(Point3D p1, Point3D p2, Point3D p3, Vector3D N)
        {
            p[0] = p1;
            p[1] = p2;
            p[2] = p3;
            n = N;
        }
    }
}
