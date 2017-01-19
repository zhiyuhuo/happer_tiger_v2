using UnityEngine;
using UnityEngine.UI;
using System.Collections;

public class TigerControl : MonoBehaviour {
	// Use this for initialization

	public GameObject angletext;
	public GameObject leg;
	private double preangle;
	private double curangle;
	private double curdir;
	private int leftLegLifts;
	private int rightLegLifts;
	private Rigidbody rb;

	private Quaternion originRotation;
	void Start () 
	{
		preangle = 0;
		curangle = 0;
		curdir = 0;
		leftLegLifts = 0;
		rightLegLifts = 0;
		originRotation = transform.rotation;
		rb = GetComponent<Rigidbody>();
	}

	// Update is called once per frame
	void Update () 
	{
		transform.rotation = originRotation;
		bool up;
		HoughControl legtrack = leg.GetComponent<HoughControl> ();
		curangle = legtrack.angle;
		curdir = legtrack.dir;
		GameObject g = GameObject.FindGameObjectWithTag ("AngleText"); 
		g.SendMessage("SetAngle", curangle);
		up = Input.GetKey (KeyCode.UpArrow);

		if (up == false) 
		{
			if (curangle >= 20 && preangle < 20) 
			{
				up = true;
				preangle = curangle;
			} 
			else 
			{
				preangle = curangle;
				up = false;
			}
		}

		if (up == true) 
		{
			//Debug.Log("upupup!");
			if (curdir == 1)
			{
				leftLegLifts ++;
				GameObject gt = GameObject.FindGameObjectWithTag ("LeftLegLifts"); 
				gt.SendMessage("SetLifts", leftLegLifts);
			}
			if (curdir == 0)
			{
				rightLegLifts ++;
				GameObject gt = GameObject.FindGameObjectWithTag ("RightLegLifts"); 
				gt.SendMessage("SetLifts", rightLegLifts);
			}
			rb.AddForce (Vector3.up * 700000 * Time.deltaTime);			
		} 
		else 
		{
			rb.AddForce (Vector3.down * 100000f * Time.deltaTime);
		}
	}
}
