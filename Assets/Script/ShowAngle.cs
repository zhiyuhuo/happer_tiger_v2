using UnityEngine;
using UnityEngine.UI;
using System.Collections;

public class ShowAngle : MonoBehaviour {

	private double angle = 0;
	// Use this for initialization
	void Start () {
		
	}
	
	// Update is called once per frame
	void Update () {
		Text t = this.GetComponents<Text> ()[0];
		t.text = "Angle: " + angle.ToString();
	}
	
	void SetAngle(double a)
	{
		angle = a;
	}
}
