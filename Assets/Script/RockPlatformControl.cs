using UnityEngine;
using System.Collections;

public class RockPlatformControl : MonoBehaviour {

	public float speed = 4;
	public bool ifGameOn;
	// Use this for initialization
	void Start () {
		ifGameOn = false;
	}
	
	// Update is called once per frame
	void Update () {
		if (ifGameOn) 
		{
			transform.Translate (Vector3.right * Time.deltaTime * speed);
		}
	}
	
	void SetState (bool game) 
	{
		ifGameOn = game;
	}
}
