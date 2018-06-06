using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Windows.Forms;

namespace kurs
{
    public partial class Form1 : Form
    {
        Bitmap bmp1;
        Bitmap bmp2;

        double dX;
        double dY;

        double kX;
        double kY;

        double K;

        private void Clear()
        {
            bmp1 = new Bitmap(pictureBox1.Width, pictureBox1.Height);
            bmp2 = new Bitmap(pictureBox2.Width, pictureBox2.Height);

            dX = 270;
            dY = 0;

            K = 15;

            Kordynaty(bmp1, K, dX, dY);
            Kordynaty(bmp2, K, dX, dY);
            pictureBox1.Image = bmp1;
            pictureBox2.Image = bmp2;
        }

        public Form1()
        {
            InitializeComponent();

            Clear();
        }

        private void button1_Click(object sender, EventArgs e)
        {
            Clear();

            double T = 1;

            double fs = 128;

            int M = 64;

            double fp = double.Parse(textBox1.Text);
            double fa = double.Parse(textBox2.Text);
            double Ap = double.Parse(textBox3.Text);
            double Aa = double.Parse(textBox4.Text);

            double f1 = double.Parse(textBox5.Text);
            double f2 = double.Parse(textBox6.Text);
            double f3 = double.Parse(textBox7.Text);
            double f4 = double.Parse(textBox8.Text);
            double f5 = double.Parse(textBox9.Text);

            double beta1 = (Math.Pow(10, 0.05 * Ap) - 1) / (Math.Pow(10, 0.05 * Ap) + 1);
            double beta2 = Math.Pow(10, -0.05 * Aa);
            double beta = Math.Min(beta1, beta2);

            double Bt = fa - fp;
            double fc = fp + Bt / 2;

            double A = -20 * Math.Log10(beta);
            double D = (A <= 21) ? (0.9222) : ((A - 7.95) / 14.36);

            double alpha = 0;
            if (A <= 21) alpha = 0;
            if (A > 21 && A <= 50) alpha = 0.5842 * Math.Pow((A - 21), 0.4) + 0.07886 * (A - 21);
            if (A > 50) alpha = 0.1102 * (A - 8.7);

            int M1 = (int)((fs * D) / Bt);

            int N1 = M1 + 1;

            double[] a = new double[M1 / 2 + 1];

            a[0] = 2 * (fc / fs);

            if (M1 % 2 == 1)
            {
                for (int i = 1; i <= M1 / 2; i++)
                {
                    a[i] = (1 / (Math.PI * (i - 0.5))) * Math.Sin(2 * Math.PI * (i - 0.5) * (fc / fs));
                }
            }
            else
            {
                for (int i = 1; i <= M1 / 2; i++)
                {
                    a[i] = (1 / (Math.PI * i)) * Math.Sin(2 * Math.PI * i * (fc / fs));
                }
            }

            Func<double, double> Io = x =>
            {
                double res = 0;

                for (int k = 1; k < 100; k++)
                {
                    res += Math.Pow(((1.0 / factorial(k)) * Math.Pow((x / 2), k)), 2);
                }
                return res + 1;
            };

            Func<double, double> wk = x => Io(alpha * (Math.Sqrt(1 - Math.Pow((2 * x) / M1, 2)))) / Io(alpha);

            Func<double, double> wb = x => 0.42 - 0.5 * Math.Cos((2 * Math.PI * x) / N1) + 0.08 * Math.Cos((4 * Math.PI * x) / N1);

            Func<double, double> wt = x => 1 - x / N1;

            Func<double, double> wh = x => alpha - (1 - alpha) * Math.Cos(2 * Math.PI * x / N1);

            double[] da = new double[M1 / 2 + 1];

            for (int i = 0; i <= M1 / 2; i++)
            {
                da[i] = a[i] * wk(i);
                //da[i] = a[i] * wb(i);
                //da[i] = a[i] * wt(i);
                //da[i] = a[i] * wh(i);
            }

            double[] h = new double[M1];

            for (int i = 0; i < M1; i++)
            {
                if (i < M1 / 2)
                {
                    h[i] = da[M1 / 2 - i];
                }

                if (i == M1 / 2)
                {
                    h[i] = da[0];
                }

                if (i > M1 / 2)
                {
                    h[i] = h[M1 - i];
                }
            }

            Func<double, double> Hw = w =>
            {
                double res = 0;

                for (int i = 1; i < M1; i++)
                {
                    res += h[i] * Math.Pow(Math.E, -i * w * M1);
                }

                return res;
            };


            textBox10.Text = "";

            for (int i = 0; i < M1; i++)
            {
                textBox10.Text += h[i].ToString("F4") + "\r\n";
            }

            Func<double, double> func1 = x =>
                Math.Cos(2 * Math.PI * f1 * x) +
                Math.Cos(2 * Math.PI * f2 * x) +
                Math.Cos(2 * Math.PI * f3 * x) +
                Math.Cos(2 * Math.PI * f4 * x) +
                Math.Cos(2 * Math.PI * f5 * x);

            FunctionDraw(func1, bmp1, 0, T, dX, dY, 35, 1, Color.Red);

            double teta = T / fs;

            int N = ((1 + (int)(T / teta)) % 2 == 0) ? (1 + (int)(T / teta)) : (1 + (int)(T / teta) + 1);

            teta = T / N;

            double[] Sd1 = new double[N];

            double[] Sd2 = new double[N];

            for (int i = 0; i < N; i++)
            {
                Sd1[i] = func1(teta * i);
            }

            double d = (Max(Sd1) - Min(Sd1)) / (M - 1);

            for (int i = 0; i < N; i++)
            {
                Sd2[i] = (int)(Sd1[i] / d);
            }

            Func<int, double> yn = k =>
            {
                double res = 0;

                for (int i = 0; ((i < h.Length) && (k - i >= 0)); i++)
                {
                    res += h[i] * Sd2[k - i];
                }

                return res;
            };

            for (int i = 0; i < N - 1; i++)
            {
                DrawLine(bmp2, new Point(teta * i, yn(i)), new Point(teta * (i + 1), yn(i + 1)), dX, dY, 1.0 / 35.0, 1.0 / d, Color.Blue);
            }

            pictureBox1.Image = bmp1;
            pictureBox2.Image = bmp2;
        }

        private double factorial(int x)
        {
            return (x == 1) ? 1 : (x * factorial(x - 1));
        }

        public void Kordynaty(Bitmap bmp, double K, double dX = 0, double dY = 0)
        {
            for (int i = (bmp.Height / 2) - (int)dY + (int)K; i < bmp.Height; i += (int)K)
            {
                for (int j = 0; j < bmp.Width; j++)
                {
                    try
                    {
                        bmp.SetPixel(j, i, Color.Silver);
                    }
                    catch { }
                }
            }

            for (int i = (bmp.Height / 2) - (int)dY - (int)K; i > 0; i -= (int)K)
            {
                for (int j = 0; j < bmp.Width; j++)
                {
                    try
                    {
                        bmp.SetPixel(j, i, Color.Silver);
                    }
                    catch { }
                }
            }

            for (int i = (bmp.Width / 2) - (int)dX + (int)K; i < bmp.Width; i += (int)K)
            {
                for (int j = 0; j < bmp.Height; j++)
                {
                    try
                    {
                        bmp.SetPixel(i, j, Color.Silver);
                    }
                    catch { }
                }
            }

            for (int i = (pictureBox1.Width / 2) - (int)dX - (int)K; i > 0; i -= (int)K)
            {
                for (int j = 0; j < bmp.Height; j++)
                {
                    try
                    {
                        bmp.SetPixel(i, j, Color.Silver);
                    }
                    catch { }
                }
            }

            for (int i = 0; i < bmp.Width; i++)
            {
                bmp.SetPixel(i, (bmp.Height / 2) - (int)dY, Color.Black);
            }

            for (int i = 0; i < bmp.Height; i++)
            {
                bmp.SetPixel((bmp.Width / 2) - (int)dX, i, Color.Black);
            }
        }

        public void FunctionDraw(Func<double, double> func, Bitmap bmp, double X1, double X2, double dX, double dY, double kX, double kY, Color color)
        {
            X1 = X1 * kX;
            X2 = X2 * kX;

            double y0 = func(X1 / kX) / kY;
            double y = 0;

            for (double i = X1; i <= X2 + 1 / K; i += 1 / K)
            {
                y = func(i / kX) / kY;

                if (i > X1)
                {
                    if (Math.Abs(Math.Abs(y) - Math.Abs(y0)) > (1 / K))
                    {
                        if (y0 > y)
                        {
                            for (double k = y; k < y0; k += (1 / K))
                                try { bmp.SetPixel((bmp.Width / 2) + (int)(i * K) - (int)dX, (bmp.Height / 2) - (int)(k * K) - (int)dY, color); }
                                catch { }
                        }
                        else
                        {
                            for (double l = y0; l < y; l += (1 / K))
                                try { bmp.SetPixel((bmp.Width / 2) + (int)(i * K) - (int)dX, (bmp.Height / 2) - (int)(l * K) - (int)dY, color); }
                                catch { }
                        }
                    }
                }

                y0 = y;

                try { bmp.SetPixel((bmp.Width / 2) + (int)(i * K) - (int)dX, (bmp.Height / 2) - (int)(y * K) - (int)dY, color); }
                catch { }
            }
        }

        double Max(double[] arr)
        {
            double max = arr[0];

            for (int i = 1; i < arr.Length; i++)
            {
                if (arr[i] > max)
                {
                    max = arr[i];
                }
            }

            return max;
        }

        double Min(double[] arr)
        {
            double min = arr[0];

            for (int i = 1; i < arr.Length; i++)
            {
                if (arr[i] < min)
                {
                    min = arr[i];
                }
            }

            return min;
        }

        private void DrawLine(Bitmap bmp, Point p1, Point p2, double dX, double dY, double kX, double kY, Color color)
        {
            double X1 = p1.X / kX;
            double X2 = p2.X / kX;

            double y0 = ((((X1 * kX - X2 * kX) * (p2.Y - p1.Y)) / (X2 * kX - X1 * kX)) + p2.Y) / kY;

            for (double i = Math.Min(X1, X2); i <= Math.Max(X1, X2); i += (1 / K))
            {
                double y = ((((i * kX - p2.X) * (p2.Y - p1.Y)) / (p2.X - p1.X)) + p2.Y) / kY;

                if (i > Math.Min(p1.X, p2.X))
                {
                    if (Math.Abs(Math.Abs(y) - Math.Abs(y0)) > (1 / K))
                    {
                        if (y0 > y)
                        {
                            for (double k = y; k < y0; k += (1 / K))
                                try { bmp.SetPixel((pictureBox1.Width / 2) + (int)(i * K) - (int)dX, (pictureBox1.Height / 2) - (int)(k * K) - (int)dY, color); }
                                catch { }
                        }
                        else
                        {
                            for (double l = y0; l < y; l += (1 / K))
                                try { bmp.SetPixel((pictureBox1.Width / 2) + (int)(i * K) - (int)dX, (pictureBox1.Height / 2) - (int)(l * K) - (int)dY, color); }
                                catch { }
                        }
                    }
                }

                y0 = y;

                try { bmp.SetPixel((bmp.Width / 2) + (int)(i * K) - (int)dX, (bmp.Height / 2) - (int)(y * K) - (int)dY, color); }
                catch { }
            }
        }
    }

    class Complex
    {
        public double Real;
        public double Imag;

        public Complex(double Re, double Im)
        {
            this.Real = Re;
            this.Imag = Im;
        }

        public double Abs()
        {
            return Math.Sqrt(Math.Pow(this.Real, 2) + Math.Pow(this.Imag, 2));
        }

        public double Fi()
        {
            return Math.Atan(this.Real / this.Imag);
        }

        public static Complex operator +(Complex A, Complex B)
        {
            return new Complex(A.Real + B.Real, A.Imag + B.Imag);
        }

        public static Complex operator *(double C, Complex A)
        {
            return new Complex(C * A.Real, C * A.Imag);
        }

        public static Complex operator *(Complex A, double C)
        {
            return new Complex(C * A.Real, C * A.Imag);
        }

        public static double Sqr(Complex c)
        {
            return (Math.Pow(c.Real, 2) + Math.Pow(c.Imag, 2));
        }

        public static double Abs(Complex c)
        {
            return Math.Sqrt(Math.Pow(c.Real, 2) + Math.Pow(c.Imag, 2));
        }

        public static double Fi(Complex c)
        {
            return Math.Atan(c.Real / c.Imag);
        }

        public override string ToString()
        {
            if (this.Imag < 0)
                return this.Real.ToString("F3") + " - j" + (-this.Imag).ToString("F3");
            return this.Real.ToString("F3") + " + j" + this.Imag.ToString("F3");
        }
    }

    public class Point
    {
        public double X;
        public double Y;

        public Point(double X, double Y)
        {
            this.X = X;
            this.Y = Y;
        }

        public override string ToString()
        {
            //return "[" + this.X.ToString() + "; " + this.Y.ToString() + "]";
            return "(" + this.X.ToString() + "; " + this.Y.ToString() + ")";
            //return "{" + this.X.ToString() + "; " + this.Y.ToString() + "}";
        }
    }
}
