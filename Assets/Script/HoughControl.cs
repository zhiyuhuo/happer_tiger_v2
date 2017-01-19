using UnityEngine;
using UnityEngine.UI;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Windows.Kinect;

public class HoughControl : MonoBehaviour {

	public GameObject showtext;
	private HoughLegControl control; 
	public double angle;
	public double dir;
	private KinectSensor sensor;
	public bool Head_at_Left;

	// Use this for initialization
	void Start () {
		control = new HoughLegControl ();
		sensor = KinectSensor.GetDefault ();
		angle = 0;
		Head_at_Left = false;
	}
	
	// Update is called once per frame
	void Update () {

		ushort[] depth = KinectManager.Instance.GetRawDepthMap ();
		CameraSpacePoint[] points = new CameraSpacePoint[512 * 424];
		sensor.CoordinateMapper.MapDepthFrameToCameraSpace (depth, points);
		control.ifHeadLeft = Head_at_Left;
		control.task (depth, points);

		angle = control.m_angle;
		dir = control.m_dir;

		//Debug.Log ("angle: " + angle.ToString ());
	}
	
}

class HoughLegControl
{
	public bool ifHeadLeft = false;
	bool ifDecideDirection = false;
	int m_dataOffset = 0;
	const int IMG_WIDTH = 512;
	const int IMG_HEIGHT = 424;
	const int IMG_SIZE = IMG_WIDTH * IMG_HEIGHT;
	
	float NEARBOUND = 1.0f;
	float FARBOUND = 3.0f;
	float LEFTBOUND = -1.2f;
	float RIGHTBOUND = 1.2f;
	
	int m_cols = 512;
	int m_rows = 424;
	
	CameraSpacePoint[] m_points = new CameraSpacePoint[IMG_SIZE];
	
	List<double> m_dirList = new List<double>(20);
	public double m_dir = 1;
	int m_count = 0;
	
	bool m_ifFloorDetected = false;
	double[] m_gd = new double[3] { 0, 0, 0 };
	int[] m_center = new int[4] { 0, 0, 0, 0 };
	
	public double m_angle = 0;
	public long m_time = 0;

	public ushort[] task(ushort[] shortDepth, CameraSpacePoint[] points)
	{
		DateTime now1 = DateTime.Now;
		
		m_points = points;
		ushort[] depth = shortDepth;
		depth = shift_arr(depth);
		m_gd = extract_floor_plane(depth);
		depth = filter_background(depth);

		/*
		if (ifHeadLeft == true) 
		{
			m_dir = 0;
		}
		else
		{
			m_dir = 1;
		}*/


		if (ifDecideDirection == false) 
		{
			m_dir = get_head_direction(depth);
			ifDecideDirection = true;
		}

        if (m_dir == 0)
        {
            depth = flip_array_horizontal(depth);
        }
		depth = resize_image(depth, 0.25);
		m_center = detect_human_center(depth);
		double angle = get_hough_angle(depth, m_center);
		Console.WriteLine("angle: " + angle.ToString());
		m_angle = angle;
		
		m_count++;
		
		DateTime now2 = DateTime.Now;
		m_time = now2.Second * 1000 + now2.Millisecond - now1.Second * 1000 - now1.Millisecond;
		
		return depth;
	}
	
	ushort[] shift_arr(ushort[] depth)
	{
		ushort[] r = new ushort[depth.Length];
		for (int i = 0; i < depth.Length; i++)
		{
			r[i] = (ushort)(depth[i] >> m_dataOffset);
		}
		
		return r;
	}
	
	double[] extract_floor_plane(ushort[] shortDepth)
	{
		ushort[] depth = shortDepth;
		List<double> ptPlane = new List<double>(shortDepth.Length);
		for (int u = 0; u < IMG_WIDTH; u += 10)
		{
			for (int v = 0; v < IMG_HEIGHT; v += 10)
			{
				if (Math.Abs(u - IMG_WIDTH / 2) < IMG_WIDTH / 8 && v > IMG_HEIGHT - IMG_HEIGHT / 8)
				{
					double[] xyz = depth_to_xyz(u, v, depth[u + v * IMG_WIDTH]);
					if (xyz[1] >= 0 && xyz[2] < 0 && xyz[0] == xyz[0])
					{
						ptPlane.Add(xyz[0]);
						ptPlane.Add(xyz[1]);
						ptPlane.Add(xyz[2]);
					}
				}
			}
		}
		double[] res = get_plane_from_point_set(ptPlane);
		m_ifFloorDetected = true;
		
		return res;
	}
	
	ushort[] filter_background(ushort[] shortDepth)
	{
		ushort[] depth = shortDepth;
		double A = m_gd[0];
		double B = m_gd[1];
		double C = m_gd[2];
		
		for (int u = 0; u < IMG_WIDTH; u++)
		{
			for (int v = 0; v < IMG_HEIGHT; v++)
			{
				int idx = v * IMG_WIDTH + u;
				double[] xyz = depth_to_xyz(u, v, depth[idx]);
				double x = xyz[0];
				double y = xyz[1];
				double z = xyz[2];
				
				if (x < LEFTBOUND
				    || x > RIGHTBOUND
				    || y < NEARBOUND
				    || y > FARBOUND
				    || z < A * x + B * y + C + 0.04
				    || z > A * x + B * y + C + 1.50)
				{
					depth[idx] = 0;
				}
				else
				{
					;
				}
			}
		}
		
		ushort[] res = depth;
		return res;
	}
	
	public double get_head_direction(ushort[] depth)
	{
		double res = 0;
		int col = IMG_WIDTH;
		int row = IMG_HEIGHT;
		double sl = 1;
		double sr = 1;
		double tl = 1;
		double tr = 1;

		for (int u = 0; u < col; u++) 
		{
			if (u < 0.45 * col)
			{
				int h = 0;
				int t = 0;
				for (int v = 20; v < row - 200; v++)
				{
					ushort p = depth[u + v * col];
					if (p > 0)
					{
						h += row - v;
						t ++;
					}
				}
				if (t > 0)
				{
					sl += h / t;
					tl ++;
				}
			}

			if (u > 0.55 * col)
			{
				int h = 0;
				int t = 0;
				for (int v = 20; v < row - 200; v++)
				{
					ushort p = depth[u + v * col];
					if (p > 0)
					{
						h += row - v;
						t ++;
					}
				}
				if (t > 0)
				{
					sr += h / t;
					tr ++;
				}
			}
		}

		sl /= tl;
		sr /= tr;

		double wl = sl;
		double wr = sr;

		if (wl > wr) 
		{
			m_dir = 0;
		}
		else
		{
			m_dir = 1;
		}

		//m_dir = wl / wr;

		return m_dir;
	}
	
	ushort[] resize_image(ushort[] depth, double p)
	{
		ushort[] res = new ushort[(int)(depth.Length * (p * p))];
		int t = 0;
		for (int i = 0; i < IMG_HEIGHT; i += (int)(1 / p))
		{
			for (int j = 0; j < IMG_WIDTH; j += (int)(1 / p))
			{
				res[t] = depth[IMG_WIDTH * i + j];
				t++;
			}
		}
		m_cols = (int)(IMG_WIDTH * p);
		m_rows = (int)(IMG_HEIGHT * p);
		
		return res;
	}
	
	int[] detect_human_center(ushort[] depth)
	{
		int[] res = new int[4] { 0, 0, 0, 0 };
		
		int ut = 0;
		int t = 0;
		
		int col = m_cols;
		int row = m_rows;
		
		for (int u = 0; u < col; u++)
		{
			for (int v = 0; v < row; v++)
			{
				if (depth[v * col + u] > 0)
				{
					ut += u;
					t++;
				}
			}
		}
		
		if (t > 0)
		{
			ut /= t;
		}
		
		
		t = 0;
		int vt = 0;
		for (int v = 0; v < row; v++)
		{
			if (depth[v * col + ut] > 0)
			{
				vt += v;
				t++;
			}
		}
		
		if (t > 0)
		{
			vt /= t;
		}
		
		
		int ve = 0;
		for (int v = row - 1; v >= 0; v--)
		{
			if (depth[v * col + ut] > 0)
			{
				ve = v;
				break;
			}
		}
		
		res[0] = ut;
		res[1] = vt;
		res[2] = ut;
		res[3] = ve;
		
		return res;
	}
	
	double get_hough_angle(ushort[] depth, int[] center)
	{
		
		double res = 0;
		int col = m_cols;
		int row = m_rows;
		int[] cen = new int[4] { 0, 0, 0, 0 };
		cen[0] = center[0];
		cen[1] = center[1];
		
		cen[0] = col - cen[0] - 1 - 5;
		cen[1] = row - cen[1] - 1;
		
		double[,] houghMap = new double[row, col];
		
		for (int v = 0; v < row; v++)
		{
			for (int u = 0; u < col; u++)
			{
				houghMap[v, u] = depth[col * (row - v - 1) + (col - u - 1)];
			}
		}
		
		double[] h = new double[90];
		
		for (int v = cen[1]; v < row; v++)
		{
			for (int u = cen[0] + 5; u < col; u++)
			{
				if (houghMap[v, u] > 0)
				{
					int angle = (int)(Math.Atan2((double)(v - cen[1]), (double)(u - cen[0])) * 180 / Math.PI);
					if (angle >= 0 && angle < 90)
					{
						double d = Math.Sqrt((v - cen[1]) * (v - cen[1]) + (u - cen[0]) * (u - cen[0]));
						d = d / ((double)(col));
						if (d < 0.1)
						{
							d = 0;
						}
						h[angle] += d;
					}
				}
			}
		}
		
		double sum = 0;
		for (int i = 0; i < h.Length; i++)
		{
			sum += h[i];
		}
		
		for (int i = 0; i < h.Length; i++)
		{
			if (h.Max() > 0)
			{
				h[i] = h[i] / sum;
			}
			else 
			{
				h[i] = 0;
			}
		}
		
		int theta = 0;
		theta = get_max_angle(h);
		res = (double)theta;
		
		return res;
	}
	
	int get_max_angle(double[] h)
	{
		int res = 0;
		
		int a = 0;
		double m = 0;
		for (int i = 0; i < h.Length; i++)
		{
			m += i * h[i];
		}
		
		a = (int)m;
		
		m = 0;
		double psum = 0;
		for (int i = a; i < h.Length; i++)
		{
			m += i * h[i];
			psum += h[i];
		}
		
		a = (int)(m / psum);
		
		res = a;
		
		return res;
	}
	
	double[] depth_to_xyz(int u, int v, ushort depth)
	{
		double[] res = new double[3] { 0, 0, 0 };
		int idx = u + v * IMG_WIDTH;
		res[0] = m_points[idx].X;
		res[2] = m_points[idx].Y;
		res[1] = m_points[idx].Z;
		
		return res;
	}
	
	double[] get_plane_from_point_set(List<double> pt)
	{
		double[] res = new double[3] { 0, 0, 0 };
		
		double xy = 0;
		double xz = 0;
		double yz = 0;
		double x2 = 0;
		double y2 = 0;
		double x = 0;
		double y = 0;
		double z = 0;
		double n = pt.Count / 3;
		for (int i = 0; i < pt.Count / 3; i++)
		{
			double xf = pt[3 * i];
			double yf = pt[3 * i + 1];
			double zf = pt[3 * i + 2];
			x += xf;
			y += yf;
			z += zf;
			
			x2 += xf * xf;
			y2 += yf * yf;
			
			xy += xf * yf;
			xz += xf * zf;
			yz += yf * zf;
		}
		
		double[,] m1Arr = new double[3, 3] {{x2, xy, x}, 
			{xy, y2, y},
			{x,  y,  n}};
		double[,] m1ArrInv = MatrixInv.InvertMatrix(m1Arr, 3);
		
		double[,] m1i = m1ArrInv;
		double[,] m2 = new double[3, 1] {{xz}, 
			{yz},
			{z}};
		
		
		res[0] = m1i[0, 0] * m2[0, 0] + m1i[0, 1] * m2[1, 0] + m1i[0, 2] * m2[2, 0];
		res[1] = m1i[1, 0] * m2[0, 0] + m1i[1, 1] * m2[1, 0] + m1i[1, 2] * m2[2, 0];
		res[2] = m1i[2, 0] * m2[0, 0] + m1i[2, 1] * m2[1, 0] + m1i[2, 2] * m2[2, 0];
		
		//Console.WriteLine("gd:" + res[0].ToString() + " " + res[1].ToString() + " " + res[2].ToString());
		return res;
	}
	
	public ushort[] flip_array_horizontal(ushort[] depth)
	{
		ushort[] res = new ushort[IMG_SIZE];
		
		for (int u = 0; u < IMG_WIDTH; u++)
		{
			for (int v = 0; v < IMG_HEIGHT; v++)
			{
				int id1 = v * IMG_WIDTH + u;
				int id2 = v * IMG_WIDTH + IMG_WIDTH - u - 1;
				res[id2] = depth[id1];
			}
		}
		
		return res;
	}
	
	byte[] get_visible_image(ushort[] depth)
	{
		byte[] res = new byte[IMG_SIZE];
		for (int i = 0; i < IMG_SIZE; ++i)
		{
			// To convert to a byte, we're mapping the depth value to the byte range.
			// Values outside the reliable depth range are mapped to 0 (black).
			res[i] = (byte)(depth[i] >= 0 ? 255 : 0);
		}
		
		return res;
	}
	
}

class MyContainer
{
	public double[][] m_double;
	public int[] m_int;
}

class MatrixInv
{
	/*
        * Perform LUP decomposition on a matrix A.
        * Return L and U as a single matrix(double[][]) and P as an array of ints.
        * We implement the code to compute LU "in place" in the matrix A.
        * In order to make some of the calculations more straight forward and to 
        * match Cormen's et al. pseudocode the matrix A should have its first row and first columns
        * to be all 0.
        * */
	public static MyContainer LUPDecomposition(double[][] A)
	{
		int n = A.Length - 1;
		/*
            * pi represents the permutation matrix.  We implement it as an array
            * whose value indicates which column the 1 would appear.  We use it to avoid 
            * dividing by zero or small numbers.
            * */
		int[] pi = new int[n + 1];
		double p = 0;
		int kp = 0;
		int pik = 0;
		int pikp = 0;
		double aki = 0;
		double akpi = 0;
		
		//Initialize the permutation matrix, will be the identity matrix
		for (int j = 0; j <= n; j++)
		{
			pi[j] = j;
		}
		
		for (int k = 0; k <= n; k++)
		{
			/*
                * In finding the permutation matrix p that avoids dividing by zero
                * we take a slightly different approach.  For numerical stability
                * We find the element with the largest 
                * absolute value of those in the current first column (column k).  If all elements in
                * the current first column are zero then the matrix is singluar and throw an
                * error.
                * */
			p = 0;
			for (int i = k; i <= n; i++)
			{
				if (Math.Abs(A[i][k]) > p)
				{
					p = Math.Abs(A[i][k]);
					kp = i;
				}
			}
			if (p == 0)
			{
				throw new Exception("singular matrix");
			}
			/*
                * These lines update the pivot array (which represents the pivot matrix)
                * by exchanging pi[k] and pi[kp].
                * */
			pik = pi[k];
			pikp = pi[kp];
			pi[k] = pikp;
			pi[kp] = pik;
			
			/*
                * Exchange rows k and kpi as determined by the pivot
                * */
			for (int i = 0; i <= n; i++)
			{
				aki = A[k][i];
				akpi = A[kp][i];
				A[k][i] = akpi;
				A[kp][i] = aki;
			}
			
			/*
                    * Compute the Schur complement
                    * */
			for (int i = k + 1; i <= n; i++)
			{
				A[i][k] = A[i][k] / A[k][k];
				for (int j = k + 1; j <= n; j++)
				{
					A[i][j] = A[i][j] - (A[i][k] * A[k][j]);
				}
			}
		}
		
		MyContainer res = new MyContainer();
		res.m_double = A;
		res.m_int = pi;
		
		return res;
	}
	
	/*
        * Given L,U,P and b solve for x.
        * Input the L and U matrices as a single matrix LU.
        * Return the solution as a double[].
        * LU will be a n+1xm+1 matrix where the first row and columns are zero.
        * This is for ease of computation and consistency with Cormen et al.
        * pseudocode.
        * The pi array represents the permutation matrix.
        * */
	public static double[] LUPSolve(double[][] LU, int[] pi, double[] b)
	{
		int n = LU.Length - 1;
		double[] x = new double[n + 1];
		double[] y = new double[n + 1];
		double suml = 0;
		double sumu = 0;
		double lij = 0;
		
		/*
            * Solve for y using formward substitution
            * */
		for (int i = 0; i <= n; i++)
		{
			suml = 0;
			for (int j = 0; j <= i - 1; j++)
			{
				/*
                    * Since we've taken L and U as a singular matrix as an input
                    * the value for L at index i and j will be 1 when i equals j, not LU[i][j], since
                    * the diagonal values are all 1 for L.
                    * */
				if (i == j)
				{
					lij = 1;
				}
				else
				{
					lij = LU[i][j];
				}
				suml = suml + (lij * y[j]);
			}
			y[i] = b[pi[i]] - suml;
		}
		//Solve for x by using back substitution
		for (int i = n; i >= 0; i--)
		{
			sumu = 0;
			for (int j = i + 1; j <= n; j++)
			{
				sumu = sumu + (LU[i][j] * x[j]);
			}
			x[i] = (y[i] - sumu) / LU[i][i];
		}
		return x;
	}
	
	/*
        * Given an nXn matrix A, solve n linear equations to find the inverse of A.
        * */
	public static double[,] InvertMatrix(double[,] Ainput, int ndim)
	{
		int n = ndim;
		double[][] A = new double[n][];
		for (int i = 0; i < n; i++)
		{
			A[i] = new double[n];
			for (int j = 0; j < n; j++)
			{
				A[i][j] = Ainput[i, j];
			}
		}
		//e will represent each column in the identity matrix
		double[] e;
		//x will hold the inverse matrix to be returned
		double[][] x = new double[n][];
		double[,] res = new double[n, n];
		for (int i = 0; i < n; i++)
		{
			x[i] = new double[A[i].Length];
		}
		/*
            * solve will contain the vector solution for the LUP decomposition as we solve
            * for each vector of x.  We will combine the solutions into the double[][] array x.
            * */
		double[] solve;
		
		//Get the LU matrix and P matrix (as an array)
		MyContainer results = LUPDecomposition(A);
		
		double[][] LU = results.m_double;
		int[] P = results.m_int;
		
		/*
            * Solve AX = e for each column ei of the identity matrix using LUP decomposition
            * */
		for (int i = 0; i < n; i++)
		{
			e = new double[A[i].Length];
			e[i] = 1;
			solve = LUPSolve(LU, P, e);
			for (int j = 0; j < solve.Length; j++)
			{
				x[j][i] = solve[j];
				res[j, i] = x[j][i];
			}
		}
		return res;
	}
}