﻿namespace STLcut
{
    partial class Form1
    {
        /// <summary>
        /// 設計工具所需的變數。
        /// </summary>
        private System.ComponentModel.IContainer components = null;

        /// <summary>
        /// 清除任何使用中的資源。
        /// </summary>
        /// <param name="disposing">如果應該處置 Managed 資源則為 true，否則為 false。</param>
        protected override void Dispose(bool disposing)
        {
            if (disposing && (components != null))
            {
                components.Dispose();
            }
            base.Dispose(disposing);
        }

        #region Windows Form 設計工具產生的程式碼

        /// <summary>
        /// 此為設計工具支援所需的方法 - 請勿使用程式碼編輯器修改
        /// 這個方法的內容。
        /// </summary>
        private void InitializeComponent()
        {
            this.openFileDialog1 = new System.Windows.Forms.OpenFileDialog();
            this.LoadSTLbutton = new System.Windows.Forms.Button();
            this.richTextBox1 = new System.Windows.Forms.RichTextBox();
            this.FindCrossTri = new System.Windows.Forms.Button();
            this.SuspendLayout();
            // 
            // openFileDialog1
            // 
            this.openFileDialog1.FileName = "openFileDialog1";
            this.openFileDialog1.Filter = "STL|*stl";
            // 
            // LoadSTLbutton
            // 
            this.LoadSTLbutton.Location = new System.Drawing.Point(12, 12);
            this.LoadSTLbutton.Name = "LoadSTLbutton";
            this.LoadSTLbutton.Size = new System.Drawing.Size(126, 23);
            this.LoadSTLbutton.TabIndex = 0;
            this.LoadSTLbutton.Text = "LoadSTL";
            this.LoadSTLbutton.UseVisualStyleBackColor = true;
            this.LoadSTLbutton.Click += new System.EventHandler(this.LoadSTLbutton_Click);
            // 
            // richTextBox1
            // 
            this.richTextBox1.Anchor = ((System.Windows.Forms.AnchorStyles)(((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Bottom) 
            | System.Windows.Forms.AnchorStyles.Right)));
            this.richTextBox1.Location = new System.Drawing.Point(331, 12);
            this.richTextBox1.Name = "richTextBox1";
            this.richTextBox1.Size = new System.Drawing.Size(608, 499);
            this.richTextBox1.TabIndex = 1;
            this.richTextBox1.Text = "";
            // 
            // FindCrossTri
            // 
            this.FindCrossTri.Location = new System.Drawing.Point(13, 67);
            this.FindCrossTri.Name = "FindCrossTri";
            this.FindCrossTri.Size = new System.Drawing.Size(125, 23);
            this.FindCrossTri.TabIndex = 2;
            this.FindCrossTri.Text = "Find cross Tri.";
            this.FindCrossTri.UseVisualStyleBackColor = true;
            this.FindCrossTri.Click += new System.EventHandler(this.FindCrossTri_Click);
            // 
            // Form1
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 12F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(951, 558);
            this.Controls.Add(this.FindCrossTri);
            this.Controls.Add(this.richTextBox1);
            this.Controls.Add(this.LoadSTLbutton);
            this.ForeColor = System.Drawing.SystemColors.ControlText;
            this.Name = "Form1";
            this.Text = "Form1";
            this.ResumeLayout(false);

        }

        #endregion

        private System.Windows.Forms.OpenFileDialog openFileDialog1;
        private System.Windows.Forms.Button LoadSTLbutton;
        private System.Windows.Forms.RichTextBox richTextBox1;
        private System.Windows.Forms.Button FindCrossTri;
    }
}

