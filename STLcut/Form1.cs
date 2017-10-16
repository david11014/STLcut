using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace STLcut
{
    public partial class Form1 : Form
    {
      
        Model M;
        List<Triangle> CrossTri;
        public Form1()
        {            
            InitializeComponent();
        }

        private void LoadSTLbutton_Click(object sender, EventArgs e)
        {
            openFileDialog1.ShowDialog();
            
            M = new Model(openFileDialog1.FileName);
            
            richTextBox1.Text = "File name: " + M.Filename.ToString() + ", Triangle count: " + M.TriNum.ToString() + "\n";
            richTextBox1.Text += "Bounding Box: " + M.BBox.MaxP.ToString() + " ~ " + M.BBox.MinP.ToString() + "\n";
        }

        private void FindCrossTri_Click(object sender, EventArgs e)
        {
            CrossTri = M.CrossTri(Convert.ToDouble(Layer.Text), 1);

            string s = "";

            foreach (Triangle T in CrossTri)
            {
                s += T.ToString() + "\n";
            }

            richTextBox1.Text += s;

        }

        private void Cut_Click(object sender, EventArgs e)
        {
            List<Point3D> crossP = new List<Point3D>();
            int rayCount = 10;
            Point3D dP = (M.BBox.MinP - M.BBox.MinP) / (rayCount + 1);
            double layer = Convert.ToDouble(Layer.Text);

            for (int i=0; i<rayCount;i++)//ray
            {
                Point3D source = new Point3D((M.BBox.MinP + dP * i).x, layer, M.BBox.MinP.z);

                foreach(Triangle T in CrossTri)
                {
                    Point3D CP = new Point3D();

                    if (IsTouch(source, new Vector3D(0, 0, 1), T, ref CP))
                        crossP.Add(CP);
                }
            }

        }

        private bool IsTouch(Point3D source, Vector3D ray, Triangle T,ref Point3D CP)
        {
            const double EPSILON = 0.0000001;
            Vector3D vertex0 = (Vector3D)T.p[0];
            Vector3D vertex1 = (Vector3D)T.p[0];
            Vector3D vertex2 = (Vector3D)T.p[0];
            Vector3D edge1, edge2, h, s, q;
            double a, f, u, v;
            edge1 = vertex1 - vertex0;
            edge2 = vertex2 - vertex0;
            h = Vector3D.Cross(ray, edge2);
            a = Vector3D.Dot(edge1, h);

            if (a > -EPSILON && a < EPSILON)
                return false;

            f = 1 / a;
            s = ray - vertex0;
            u = f * (Vector3D.Dot(s,h));
            if (u < 0.0 || u > 1.0)
                return false;

            q = Vector3D.Cross(s,edge1);
            v = f * Vector3D.Dot(ray, q);
            if (v < 0.0 || u + v > 1.0)
                return false;

            // At this stage we can compute t to find out where the intersection point is on the line.
            double t = f * Vector3D.Dot(edge2, q);
            if (t > EPSILON) // ray intersection
            {
                CP = (Point3D)(ray + (ray.unit() * (t * ray.length())));
                return true;
            }
            else // This means that there is a line intersection but not a ray intersection.
                return false;

        }
    }

    
    
   
}
