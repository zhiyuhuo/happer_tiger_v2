using UnityEngine;
using UnityEngine.UI;
using System.Collections;

public class ChangeScore : MonoBehaviour {

	private int score = 0;
	// Use this for initialization
	void Start () {
	
	}
	
	// Update is called once per frame
	void Update () {
		Text t = this.GetComponents<Text> ()[0];
		t.text = "Score: " + score.ToString();
	}

	void AddScore()
	{
		score++;
	}
}
