using UnityEngine;
using UnityEngine.UI;
using System.Collections;

public class ShowLeftLegLifts : MonoBehaviour {

	private double lifts;
	// Use this for initialization
	void Start () {
		lifts = 0;
	}
	
	// Update is called once per frame
	void Update () {
		Text t = this.GetComponents<Text> ()[0];
		t.text = "Left Leg Lifts: " + lifts.ToString();
	}

	void SetLifts(int a)
	{
		lifts = a;
	}
}
