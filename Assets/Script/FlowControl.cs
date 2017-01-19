using UnityEngine;
using System.Collections;

public class FlowControl : MonoBehaviour {
	
	public GameObject rockplatform;
	public GameObject coins;
	public GameObject columns;
	public GameObject background;
	public GameObject FirstBG;
	
	public float speed = 4;
	private Quaternion originBackgroundrotation;
	private Vector3 PrevRPpose1;
	private Vector3 PrevRPpose2;

	public bool ifTriggerGameStart;
	public bool ifGameStart;

	// Use this for initialization
	void Start () 
	{
		ifTriggerGameStart = true;
		ifGameStart = false;
	}
	
	// Update is called once per frame
	void Update () 
	{
		Debug.Log (ifTriggerGameStart.ToString ());
		if (ifTriggerGameStart == true) 
		{
			originBackgroundrotation = FirstBG.transform.rotation;
			PrevRPpose1 = rockplatform.transform.position;
			PrevRPpose2 = rockplatform.transform.position;
			ifTriggerGameStart = false;
			ifGameStart = true;
		}
		Debug.Log (ifGameStart.ToString ());
		if (ifGameStart == true) 
		{
			GameObject[] coinsarchset = GameObject.FindGameObjectsWithTag ("CoinsArch"); 
			foreach (GameObject c in coinsarchset) {
				c.SendMessage ("SetState", true);
			}
			GameObject[] columnsarrayset = GameObject.FindGameObjectsWithTag ("ColumnArray"); 
			foreach (GameObject c in columnsarrayset) {
				c.SendMessage ("SetState", true);
			}
			GameObject[] backgroundset = GameObject.FindGameObjectsWithTag ("Background"); 
			foreach (GameObject b in backgroundset) {
				b.SendMessage ("SetState", true);
			}
			GameObject[] rockplatformset = GameObject.FindGameObjectsWithTag ("RockPlatform"); 
			foreach (GameObject r in rockplatformset) {
				r.SendMessage ("SetState", true);
			}

			if (Input.GetKey ("escape")) {
				Application.Quit ();
			}
			float dx1 = PrevRPpose1.x - rockplatform.transform.position.x;
			if (dx1 < -8.0f) {
				PrevRPpose1 = rockplatform.transform.position;
				Instantiate (coins, new Vector3 (-14f, 2, 0), Quaternion.identity);
				Instantiate (columns, new Vector3 (-20f, 0, 0), Quaternion.identity);
			}

			float dx2 = PrevRPpose2.x - rockplatform.transform.position.x;
			//Debug.Log ("dx2: " + dx2);
			if (dx2 < - 50) {
				PrevRPpose2 = rockplatform.transform.position;
				Instantiate (background, new Vector3 (-50, 5, -6), originBackgroundrotation);
			}
		}
	}
}
