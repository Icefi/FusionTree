#ifndef FUSIONTREE_HPP
#define FUSIONTREE_HPP

#define None -1

class FusionTree
{  
	public:
		/* FusionTree fields */

		unsigned int w;
		unsigned int r;

		/* FusionTree Node */
		class Node
		{
			public:
				/* Fusion Tree Node fields */
				bool isLeaf;

				long int* keys; // key_count
				Node** children; // key_count + 1 by B-tree properties
				unsigned int key_count; // number of keys

				long int m; // mysterious constant

				long int bs; // distinguihsing bits
				long int ms; // m_j bits

				long int sketches;
				long int sketch_maskl; //mask q
				long int sketch_maskh; //mask b

				long int b_mask;
				long int bm_mask;
					
				unsigned int sketch_gap;

				unsigned int r; // the maximum number of keys

				/* Fusion Tree Node methods */
				Node (unsigned int);

				~Node ();

				void precompute (unsigned int);

				void info ();

				void merge (long int);

				void fill (long int);
		};

	Node* root;

	/* FusionTree functions */
	
	FusionTree ();

	~FusionTree ();
		
	void splitChild (Node*, int);

	void traverse (Node*);

	void traverse ();

	void init (Node*);

	void init ();

	void fill (Node*, long int);

	void removeFromLeaf (Node*, long int);

	void removeFromNonLeaf (Node*, long int);

	void remove (long int, Node*);

	void precompute (Node*);

	int parallelComparison (Node*, long int);

	long int getDiffBits (Node*);

	long int sketchApprox (Node*, long int);

	long int successor (Node*, long int);

	long int successor (long int);

	long int predecessor (Node*, long int);

	long int predecessor (long int);

	long int get_m_bits (long int);

	Node* lookUp (long int, Node*);

	Node* insertNormal (Node*, long int);

	Node* insert (long int);

	/* Static Fusion Tree Functions */

	static long int get_m (long int);

	static long int getMask (long int);

	static long int getComboMask (long int, long int);
}; 

#endif
