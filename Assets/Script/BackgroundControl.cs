using UnityEngine;
using System.Collections;

public class BackgroundControl : MonoBehaviour {
	
	public bool ifGameOn;
	public float speed = 4;
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
		if (this.transform.position.x > 50) 
		{
			Destroy (this.gameObject);
		}	
	}

	void SetState (bool game) 
	{
		ifGameOn = game;
	}
}
