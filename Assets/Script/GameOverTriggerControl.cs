using UnityEngine;
using System.Collections;

public class GameOverTriggerControl : MonoBehaviour {
	

	// Use this for initialization
	void Start () {
	
	}

	void OnCollisionEnter(Collision other) 
	{
		Debug.Log ("Collision ENTERED!");
		if (other.gameObject.name == "Tiger") 
		{
			Application.Quit();
			Debug.LogError("GAME OVER!!!");
		} 

			string[] movabletags = new string[]{"CoinsArch", "ColumnArray", "RockPlatform", "Background"};

		foreach (string tagstr in movabletags) 
		{
			GameObject[] gos = GameObject.FindGameObjectsWithTag (tagstr); 
			foreach (GameObject g in gos) 
			{
				g.SendMessage ("SetState", false);
			}
		}
	
	}

	void Update()
	{

	}
	
}
