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
        string Filename;

        Model M;
        public Form1()
        {
            Filename = "";
            M = new Model();
            InitializeComponent();
        }

        private void LoadSTLbutton_Click(object sender, EventArgs e)
        {
            openFileDialog1.ShowDialog();
            Filename = openFileDialog1.FileName;
            M.Filename = Filename;
            M.LoadSTLBinnary();
        }

    }

    
   
}
