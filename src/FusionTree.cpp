#include <iostream>
#include <cmath>

#include "FusionTree.hpp"

/* Fusion Tree Node Functions */

FusionTree::Node::Node (unsigned int max_keys = 0)
{
	key_count = 0;
	isLeaf = true;
	r = max_keys;

	keys = new long int [max_keys];
	children = new Node* [max_keys + 1];
}

FusionTree::Node::~Node ()
{
	delete keys;

	for (unsigned int i = 0; i <= key_count; i++)
		delete children [i];
}

void FusionTree::Node::info ()
{
	std::cout << "    | Node: " << this << "\n";
	std::cout << "    | bs = " << bs << "\n";
	std::cout << "    | m = " << m << "\n";
	std::cout << "    | sketch = " << sketches << "\n";
	std::cout << "    | Keys: ";
	

	for (unsigned int i = 0; i < key_count; i++)
		std::cout << keys[i] << " ";

	std::cout << "\n";
}

void FusionTree::Node::merge (long int ind)
{
	Node* child = children[ind];
	Node* sibling = children [ind + 1];

	child->keys[r / 2 - 1] = keys[ind];

	for (unsigned int i = 0; i < sibling->key_count; i++)
		child->keys[i + r / 2] = sibling->keys[i];

	if (!child->isLeaf) {
		for (unsigned int i = 0; i <= sibling->key_count; i++)
			child->children[i + r / 2] = sibling->children[i];

		for (unsigned int i = ind + 1; i < key_count; i++)
			keys[i - 1] = keys[i];

		for (unsigned int i = ind + 2; i <= key_count; i++)
			children[i - 1] = children [i];

		child->key_count += sibling->key_count + 1;
		key_count--;

		delete sibling;
	}
}

/* Fusion Tree Functions */

FusionTree::FusionTree () : w (243), r (3) { this->root = new Node (r); }

FusionTree::~FusionTree () { delete root; }

void FusionTree::precompute (Node* node)
{
	if (node->key_count != 0) {
		
		node->bs = getDiffBits (node);
		node->ms = get_m_bits (node->bs);
		
		node->m = FusionTree::get_m (node->ms);

		node->b_mask = getMask (node->bs);
		node->bm_mask = getComboMask (node->bs, node->m);

		node->sketches = 0;
		node->sketch_gap = node->key_count * node->key_count * node->key_count + 1L;

		node->sketch_maskh = 0;
		node->sketch_maskl = 0;

		long int sketch;

		for (unsigned int i = 0; i < node->key_count; i++) {
			sketch = sketchApprox(node, node->keys[i]);
			node->sketches |= (sketch | (1L << node->sketch_gap)) << i * (node->sketch_gap + 1L);
			node->sketch_maskl |= 1L << i * (node->sketch_gap + 1L);
			node->sketch_maskh |= (1L << node->sketch_gap) << i * (node->sketch_gap + 1L);
		}
	}
}

void FusionTree::traverse (Node* node)
{
	if (node == nullptr) return;

	node->info ();

   	std::cout << "\n";
   	for (unsigned int i = 0; i <= node->key_count; i++)
		traverse (node->children[i]);
}

void FusionTree::traverse ()
{
	std::cout << "FUSION TREE INFO:\n\n";
	std::cout << "    |R is " << this->r << std::endl;
	std::cout << "    |W if " << this->w << std::endl;
	std::cout << "\nFUSION TREE NODES:\n\n";
	//lookUp (24);
	traverse (root);
}

void FusionTree::splitChild (Node* node, int x)
{
	Node* z = new Node (r);
	Node* y = node->children[x];

	unsigned int pos_key = (r / 2);

	z->key_count = r - pos_key - 1;

	for (unsigned int i = 0; i < z->key_count; i++) {
		z->keys[i] = y->keys[pos_key + i + 1];
		y->keys[pos_key + i + 1] = None;
	}

	if (!y->isLeaf) {
		for (unsigned int i = 0; i <= z->key_count; i++)
			z->children[i] = y->children[pos_key + i + 1];
	}
			
	y->key_count = r - z->key_count - 1;

	node->keys[x] = y->keys[pos_key];
			
	y->keys[pos_key] = None;

	node->children[x + 1] = z;
	node->key_count++;
}

void FusionTree::init (Node* node)
{
	precompute (node);
	
	if (!node->isLeaf) {
		for (unsigned int i = 0; i <= node->key_count; ++i)
			init (node->children[i]);
	}
}

void FusionTree::init ()
{
	init (root);
}

int FusionTree::parallelComparison (Node* node, long int key)
{
	long int sketch = sketchApprox (node, key);

	long int sketch_long = sketch * node->sketch_maskl;

	long int res = node->sketches - sketch_long;

	res &= node->sketch_maskh;

	int i = 0;
	while ((1 << i) < res) i++;

	i++;
	node->sketch_gap = pow (node->key_count, 3) + 1;

	return node->key_count - (i / node->sketch_gap);
}

void FusionTree::fill (Node* node, long int ind)
{
	if ((ind != 0) && (node->children[ind - 1]->key_count >= (r /2))) {
		predecessor (node->children[ind - 1], ind);
	} else if ((ind != node->key_count) && (node->children[ind + 1]->key_count >= (r / 2))) {
		successor (node->children[ind + 1], ind);
	} else {
		if (ind != node->key_count) node->merge (ind);
		else node->merge (ind - 1);
	}
}

void FusionTree::removeFromLeaf (Node* node, long int ind)
{
	for (unsigned int i = ind + 1; i < node->key_count; ++i)
		node->keys[i - 1] = node->keys [i];
	
	node->key_count--;
}

void FusionTree::removeFromNonLeaf (Node* node, long int ind)
{
	long int k = node->keys[ind];

	if (node->children[ind]->key_count >= (r / 2)) {

		long int pred = predecessor (node->children[ind], k);
		remove (pred, node->children[ind]);
	} else if (node->children[ind + 1]->key_count >= (r / 2)) {

		long int succ = successor (node->children[ind + 1], k);
		remove (succ, node->children[ind + 1]);
	} else {
		node->merge (ind);
		remove (k, node->children[ind]);
	}
}

void FusionTree::remove (long int key, Node* node = nullptr) 
{
	Node* search = lookUp (key, node);

	if (search == nullptr) return;

	long int pos = parallelComparison (search, key);

	if (search->isLeaf) { removeFromLeaf (search, pos); }
	else { removeFromNonLeaf (search, pos); }

	precompute (search);
}

long int FusionTree::getDiffBits (Node* node)
{
	long int b_bits = 0;
	long int res = 0;

	unsigned int temp;

	for (unsigned int i = 0; i < node->key_count; i++) {
		if (node->keys[i] == None) break;
		for (unsigned int j = 0; j < i; j++) {
			temp = w;
			
			while (((node->keys[i] & (1L << temp)) == (node->keys[j] & (1L << temp))) && (temp >= 0L)) temp--;
				
			if (temp >= 0L) b_bits |= 1L << temp;
		}
	}

	for (unsigned int i = 0; i < 64; i++) {
		if ((b_bits & (1L << i)) != 0)
			res |= i;
	}

	return b_bits;
}

long int FusionTree::get_m_bits (long int b_bits)
{
	long int t, i, j, k;
	long int m_bits = 0;
	long int mt;

	bool flag;
	
	for (t = 0; t < 64; t++) {
		
		mt = 0;
		flag = true;

		while (flag) {
			flag = false;

			for (i = 0; i < 64; i++) {
				if (flag) break;

				for (j = 0; i < 64; i++) {
					if (flag) break;

					for (k = 0; k < t; k++) {
						if (mt == ((b_bits & (1L << i)) - (b_bits & (1L << j)) + (m_bits & (1L << k)))) {
							flag = true;
							break;
						}
						//std::cout << mt << "\n";
					}
				}
			}

			if (flag == true) mt++;
		}

		m_bits |= (1L << mt);
	}

	return m_bits;
}

long int FusionTree::sketchApprox (Node* node, long int x)
{
	long int xx = x & node->b_mask;
	long int res = xx * node->m;
	res &= node->bm_mask;
	return res;
}

long int FusionTree::successor (Node* node, long int key)
{
	long int res;
	unsigned int pos;

	if (node == nullptr) node = root;

	if (node->key_count == 0)
		return (node->isLeaf) ? -1 : successor (node->children[0], key);
	
	if (node->keys[0] >= key) {
		if (!node->isLeaf) {
			res = successor (node->children [0], key);
			if (res == -1)
				return node->keys[0];
			else
				return (node->keys[0] <= res) ? node->keys[0] : res;
		} else {
			return node->keys[0];
		}
	}

	if (node->keys[node->key_count - 1] < key)
		return (node->isLeaf) ? -1 : successor (node->children[node->key_count], key);
	
	pos = parallelComparison (node, key);
	//std::cout << "pos = " << pos << std::endl;
	
	if (pos == 0) pos = 1;

	long int x = (node->keys[pos - 1] > node->keys[pos]) ? node->keys[pos - 1] : node->keys[pos];
	long int common_prefix = 0;
	int i = w;

	while ((i >= 0) && ((x & (1 << i)) == (key & (1 << i)))) {
		common_prefix |= x & (1 << i);
		i--;
	}

	if (i == -1) return x;

	long int temp = common_prefix | (1 << i);
	
	pos = parallelComparison (node, temp);

	if (node->isLeaf) return node->keys[pos];
	else {
		res = successor (node->children[pos], key);
		if (res == -1)
			return node->keys[pos];
		else
			return res;
	}
}

long int FusionTree::successor (long int key) { return successor (root, key); }

long int FusionTree::predecessor (Node* node, long int key)
{
	long int res;
	unsigned int pos;

	if (node == nullptr) node = root;

	if (node->key_count == 0)
		return (node->isLeaf) ? -1 : predecessor (node->children[0], key);
	
	if (node->keys[0] > key) {
		if (!node->isLeaf) return predecessor (node->children[0], key);
		else return -1;
	}

	if (node->keys[node->key_count - 1] <= key) {
		if (node->isLeaf) return node->keys[node->key_count - 1];
		else {
			res = predecessor (node->children[node->key_count], key);
			return res > node->keys[node->key_count - 1] ? res : node->keys[node->key_count - 1];
		}
	}

	pos = parallelComparison (node, key);

	if (pos == 0) pos = 1;

	long int x = node->keys[pos];
	long int common_prefix = 0;

	int i = w;

	while ((i >= 0) && ((x & (1 << i)) == (key & (1 << i)))) {
		common_prefix |= x & (1 << i);
		i--;
	}

	if (i == -1) return x;

	long int temp = common_prefix | ((1 << i) - 1);
	pos = parallelComparison (node, temp);

	if (pos == 0) {
		if (node->isLeaf) return node->keys[pos];

		res = predecessor (node->children[1], key);
		if (res == -1) return node->keys[pos];
		else return res;
	}

	if (node->isLeaf) return node->keys[pos - 1];
	else {
		res = predecessor (node->children[pos], key);
		if (res == -1) return node->keys[pos - 1];
		else return res;
	}
}

long int FusionTree::predecessor (long int key) { return predecessor (root, key); }

FusionTree::Node* FusionTree::lookUp (long int key, Node* node = nullptr)
{
	node = (node == nullptr) ? root : node;
	
	long int pos = parallelComparison (node, key);

	if (node->keys[pos] == key) return node;
	else if (!node->isLeaf) return lookUp (key, node->children[pos]);
	else return nullptr;
}

FusionTree::Node* FusionTree::insertNormal (Node* node, long int key)
{
	int i = node->key_count;

	if (node->isLeaf) {

		while ((i >= 1) && (key < node->keys[i - 1])) {
			node->keys[i] = node->keys[i - 1];
			i--;
		}

		node->keys[i] = key;
		node->key_count++;
		return node;

	} else {

		while ((i >= 1) && (key < node->keys[i - 1])) i--;

		if (node->children[i]->key_count == r) {
			splitChild (node, i);

			if (key > node->keys[i]) i++;
		}

		return insertNormal (node->children[i], key);
	}
}

FusionTree::Node* FusionTree::insert (long int key)
{
	if (root->key_count == r) {

		Node* temp_node = new Node (r);
		temp_node->isLeaf = 0;
		temp_node->key_count = 0;
		temp_node->children[0] = root;

		root = temp_node;
		splitChild (temp_node, 0);
		return insertNormal (temp_node, key);

	} else {
		return insertNormal (root, key);
	}
}

/* Static Fusion Tree functions */
long int FusionTree::get_m (long int m_bits)
{
	long int m = 0;
	
    for (int i = 0; i < 64; i++) {
    	m |= (m_bits & (1L << i)) ? i : 0;
	}

	std::cout << "m_j = " << m_bits;
	std::cout << "\nm = " << m << std::endl;
	return m;
}

long int FusionTree::getMask (long int x)
{
	long int res = 0;

	for (int i = 0; i < 64; i++)
		res |= x & (1L << i);

	//std::cout << "Res is " << res << "\nX is " << x << std::endl;
	return res;
}

long int FusionTree::getComboMask (long int x, long int y) { return getMask (x + y); }


