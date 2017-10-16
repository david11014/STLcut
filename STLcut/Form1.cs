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
            List<Triangle> CrossTri = M.CrossTri(10, 1);

            string s = "";

            foreach (Triangle T in CrossTri)
            {
                s += T.ToString() + "\n";
            }

            richTextBox1.Text += s;

        }
    }

    
   
}
