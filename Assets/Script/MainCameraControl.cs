using UnityEngine;
using System.Collections;

public class MainCameraControl : MonoBehaviour {

	public GameObject tiger;
	// Use this for initialization
	void Start () {

	}
	
	// Update is called once per frame
	void Update () {
		Vector3 pos = new Vector3 (tiger.transform.position.x, transform.position.y, transform.position.z);
		transform.position = pos;
	}
}
