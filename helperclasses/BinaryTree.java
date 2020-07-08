package helperclasses;

import java.util.ArrayList;
import java.util.Random;

//This code is adapted from https://stackoverflow.com/questions/828398/how-to-create-a-binary-tree and https://algorithms.tutorialhorizon.com/print-a-path-from-root-to-node-in-binary-tree/
public class BinaryTree {
	
	private Node datum;
    private BinaryTree left;
    private BinaryTree right;
    
    
	public BinaryTree(Node[] values) {this(values, 0);}

	private BinaryTree(Node[] values, int index)
	{
	   Load(this, values, index);
	}
	
	//generates a binary tree with N nodes containing 1, 2, ..., N
	public static BinaryTree generateTree(int N) {
		Node[] values = new Node[N];
		for (int i = 0; i < N; i++) {
			values[i] = new Node(i + 1);
		}
		return new BinaryTree(values);
	}

	 private void Load(BinaryTree tree, Node[] values, int index)
	 {
	       this.datum = values[index];
	       if (index * 2 + 1 < values.length)
	       {
	           this.left = new BinaryTree(values, index * 2 + 1);
	       }
	       if (index * 2 + 2 < values.length)
	       {
	           this.right = new BinaryTree(values, index * 2 + 2);
	       }
	 }  
	  
	 //gets the path to a certain node -- O(N). Note that the path is reversed.
	 public ArrayList<Node> getPath(Node dest) {
		 ArrayList<Node> path = new ArrayList<Node>();
		 getPathHelper(path, this, dest);
		 return path;
	 }
	 
	 
	 public Node getRandomLeafNode() {
		 return getRandomLeafNodeHelper(this);
	 }
	 
	 private Node getRandomLeafNodeHelper(BinaryTree root) {	 
		 //random: 0 = left, 1 = right
		 if (root.left == null && root.right == null) return root.datum;
		 else if (root.left == null)
			 return getRandomLeafNodeHelper(root.right);
		 else if (root.right == null)
			 return getRandomLeafNodeHelper(root.left);
		 else {
			 int random = new Random().nextInt(2);
			 if (random == 0)
				 return getRandomLeafNodeHelper(root.left);
			 else
				 return getRandomLeafNodeHelper(root.right);
		 } 
	 }
	  
	 private boolean getPathHelper(ArrayList<Node> path, BinaryTree root, Node dest) {
		  
		 if (root == null) return false;
		 if(root.datum.equals(dest) ||getPathHelper(path, root.left,dest) || getPathHelper(path, root.right,dest)){
			//System.out.print("  " + root.data);
			path.add(root.datum);
			return true;
		}
		
		return false;  
	  }
	 
	 //This was found on https://stackoverflow.com/questions/4965335/how-to-print-binary-tree-diagram (to be used for debugging purposes)
	 public static void printBinaryTree(BinaryTree root, int level){
		    if(root==null)
		         return;
		    printBinaryTree(root.right, level+1);
		    if(level!=0){
		        for(int i=0;i<level-1;i++)
		            System.out.print("|\t");
		            System.out.println("|-------"+root.datum);
		    }
		    else
		        System.out.println(root.datum);
		    printBinaryTree(root.left, level+1);
		} 
	}
