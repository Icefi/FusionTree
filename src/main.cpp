#include <iostream>
#include "FusionTree.hpp"

int main ()
{
	FusionTree tree;

	tree.insert (1);
	tree.insert (5);
	tree.insert (15);
	tree.insert (16);
	tree.insert (20);
	tree.insert (25);
	tree.insert (4);

	tree.init ();

	tree.traverse ();

	std::cout << "21 -> " << tree.successor (21) << "\n";

	return 0;
} 
