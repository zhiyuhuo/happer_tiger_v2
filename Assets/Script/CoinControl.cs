using UnityEngine;
using System.Collections;

public class CoinControl : MonoBehaviour {
	
	public AudioClip clip;
	// Use this for initialization
	void Start () {
	
	}
	
	void OnTriggerEnter(Collider other) 
	{
		//Debug.Log ("TRIGGER ENTERED!");
		if (other.gameObject.tag == "Player") 
		{
			GameObject g = GameObject.FindGameObjectWithTag ("ScoreText"); 
			g.SendMessage("AddScore");
			//Debug.Log ("COLLECTED A COIN! (In Coin Controller)");
			GetComponent<AudioSource> ().PlayOneShot (clip);
			AudioSource.PlayClipAtPoint (clip, transform.position);
			Destroy (this.gameObject);
		} 
	}
	
	// Update is called once per frame
	void Update () {
		transform.Rotate (Vector3.forward, 100.0f*Time.deltaTime);
		//Debug.Log (this.transform.parent.transform.position.x);
	}
}
