#pragma once

template<class elementClass, class keyClass, class nodeClass> class AVLNodeLogic
{
private:
    static void		 initialize(nodeClass self);

    static void		 insert(
					nodeClass							self,
				    AVLTree<elementClass,keyClass>		*intoTree,
				    elementClass						*element);

    static void		 remove(
					nodeClass							self,
				    AVLTree<elementClass,keyClass>		*fromTree);

    static unsigned		 checkAndReturnDepth(
					nodeClass							self,
				    unsigned							*countedElements);

    static int			 inTree(nodeClass self);    

	static keyClass	 getKey(nodeClass self)
					{
						return self->getElement()->getKey();
					}

    static nodeClass
				findFirstLessThan(
					nodeClass							self,
				    keyClass							 key);

    static nodeClass
				findFirstLessThanOrEqualTo(
					nodeClass							self,
				    keyClass							 key);

    static nodeClass
				findFirstGreaterThan(
					nodeClass							self,
				    keyClass							 key);

    static nodeClass
				findFirstGreaterThanOrEqualTo(
					nodeClass							self,
				    keyClass							 key);

    static void		 rightAdded(
					nodeClass							self,
				    AVLTree<elementClass,keyClass>		*tree);

    static void		 leftAdded(
					nodeClass							self,
				    AVLTree<elementClass,keyClass>		*tree);

    static void		 singleRotate(
					nodeClass							self,
				    AVLTree<elementClass,keyClass>		*tree,
				    nodeClass							child,
				    AVLBalance			 				 whichSide);

    static void		 doubleRotate(
					nodeClass							self,
				    AVLTree<elementClass,keyClass>		*tree,
				    nodeClass							child,
				    nodeClass							grandchild,
				    AVLBalance			 				 whichSide);

    static void		 gotOneShorter(
					nodeClass							self,
				    AVLTree<elementClass,keyClass>		*tree,
				    AVLBalance			 				 whichSide);

    friend class AVLTree<elementClass,keyClass>;
};

template<class elementClass, class keyClass, class nodeClass> void
AVLNodeLogic<elementClass,keyClass,nodeClass>::initialize(
	nodeClass				self)
{
    self->setBalance(AVLNew);
	self->setLeft(NULL);
	self->setRight(NULL);
	self->setParent(NULL);
    self->setElement(NULL);
}


//****************************************************************************
//*                                                                          *
//* Function:  findFirstLessThan                                             *
//*                                                                          *
//* Syntax:    AVLElement * findFirstLessThan(                               *
//*                         elementClass * element)                          *
//*                                                                          *
//* Input:     elementClass * element:                                       *
//*              A pointer to an element to compare against while searching. *
//*                                                                          *
//* Output:    AVLElement *:                                                 *
//*              The element in the tree that has a value less than or equal *
//*              to the one specified, or NULL on failure.                   *
//*                                                                          *
//* Synopsis:  This function finds the element in the tree that has a value  *
//*            less than or equal to the one specified.                      *
//*                                                                          *
//**************************************************************************** 
template<class elementClass, class keyClass, class nodeClass> nodeClass
AVLNodeLogic<elementClass,keyClass,nodeClass>::findFirstLessThan(
					nodeClass				self,
					keyClass 				key)
{
    nodeClass retVal = NULL;
    
    if (self->getElement()->getKey() < key) {
		// The current element is smaller than the one specified.
		// This might be it, but try to find a bigger one.

		if (self->getRight() != NULL) {
		    retVal = findFirstLessThan(self->getRight(), key);
		}
	
		// If nothing below us (to the right) was found, then we are the
		// next smallest one.
		if (retVal == NULL) {
		    return self;
		} else {
		    return retVal;
		}
    } else {
		// The current element is strictly bigger than the one specified.
		// We have to find a smaller one.
		if (self->getLeft() != NULL) {
		    return findFirstLessThan(self->getLeft(), key);
		} else {
		    return NULL;
		}
    }
}
//****************************************************************************
//*                                                                          *
//* Function:  findFirstLessThanOrEqualTo                                    *
//*                                                                          *
//* Syntax:    AVLElement * findFirstLessThanOrEqualTo(                      *
//*                         elementClass * element)                          *
//*                                                                          *
//* Input:     elementClass * element:                                       *
//*              A pointer to an element to compare against while searching. *
//*                                                                          *
//* Output:    AVLElement *:                                                 *
//*              The element in the tree that has a value less than or equal *
//*              to the one specified, or NULL on failure.                   *
//*                                                                          *
//* Synopsis:  This function finds the element in the tree that has a value  *
//*            less than or equal to the one specified.                      *
//*                                                                          *
//**************************************************************************** 
template<class elementClass, class keyClass, class nodeClass> nodeClass
AVLNodeLogic<elementClass,keyClass,nodeClass>::findFirstLessThanOrEqualTo(
					nodeClass				self,
					keyClass 				key)
{
    nodeClass retVal = NULL;
    
    if (self->getElement()->getKey() == key) {
		// We have a direct match (equal to).  It takes precidence over the
		// "first less than" part.
		return self;
    }

    if (self->getElement()->getKey() < key) {
		// The current element is smaller than the one specified.
		// This might be it, but try to find a bigger one.

		if (self->getRight() != NULL) {
		    retVal = findFirstLessThanOrEqualTo(self->getRight(), key);
		}
	
		// If nothing below us (to the right) was found, then we are the
		// next smallest one.
		if (retVal == NULL) {
		    return self;
		} else {
		    return retVal;
		}
    } else {
		// The current element is bigger than the one specified.
		// We have to find a smaller one.
		if (self->getLeft() != NULL) {
		    return findFirstLessThanOrEqualTo(self->getLeft(), key);
		} else {
		    return NULL;
		}
    }
}

//****************************************************************************
//*                                                                          *
//* Function:  findFirstGreaterThan                                          *
//*                                                                          *
//* Syntax:    AVLElement * findFirstGreaterThan(elementClass * element)     *
//*                                                                          *
//* Input:     elementClass * element:                                       *
//*              A pointer to an element to compare against while searching. *
//*                                                                          *
//* Output:    AVLElement *:                                                 *
//*              The element in the tree that has a vlaue greater than the   *
//*              one specified, or NULL on failure.                          *
//*                                                                          *
//* Synopsis:  This function finds the element in the tree that has a value  *
//*            greater than the one specified.                               *
//*                                                                          *
//**************************************************************************** 
template<class elementClass, class keyClass, class nodeClass> nodeClass
AVLNodeLogic<elementClass,keyClass,nodeClass>::findFirstGreaterThan(
			nodeClass				self,
			keyClass 				key)
{
    AVLElement<elementClass,keyClass> * retVal = NULL;
    
    if (self->getElement()->getKey() > key) {
		// The current element is bigger than the one specified.
		// This might be it, but try to find a smaller one.
		if (self->getLeft() != NULL) {
		    retVal = findFirstGreaterThan(self->getLeft(), key);
		}

		// If nothing below us (to the left) was found, then we are the
		// next biggest one.
		if (retVal == NULL) {
		    return self;
		} else {
		    return retVal;
		}
    } else {
		// The current element is smaller than (or equal) the one specified.
		// We have to find a bigger one.
		if (self->getRight() != NULL) {
	    	return findFirstGreaterThan(self->getRight(), key);
		} else {
		    return NULL;
		}
    }
}

//****************************************************************************
//*                                                                          *
//* Function:  findFirstGreaterThanOrEqualTo                                 *
//*                                                                          *
//* Syntax:    AVLElement * findFirstGreaterThanOrEqualTo(elementClass * element)
//*                                                                          *
//* Input:     elementClass * element:                                       *
//*              A pointer to an element to compare against while searching. *
//*                                                                          *
//* Output:    AVLElement *:                                                 *
//*              The element in the tree that has a vlaue greater than or    *
//*              equal to the one specified, or NULL on failure.             *
//*                                                                          *
//* Synopsis:  This function finds the element in the tree that has a value  *
//*            greater than or equal to the one specified.                   *
//*                                                                          *
//**************************************************************************** 
template<class elementClass, class keyClass, class nodeClass> nodeClass
AVLNodeLogic<elementClass,keyClass,nodeClass>::findFirstGreaterThanOrEqualTo(
			nodeClass				self,
			keyClass 				key)
{
    if (self->getElement()->getKey() == key) {
		// We have a direct match (equal to).  It takes precidence over the
		// "first less than" part.
		return self;
    }

    AVLElement<elementClass,keyClass> *retVal = NULL;
    
    if (self->getElement()->getKey() > key) {
		// The current element is bigger than the one specified.
		// This might be it, but try to find a smaller one.
		if (self->getLeft() != NULL) {
		    retVal = findFirstGreaterThanOrEqualTo(self->getLeft(), key);
		}

		// If nothing below us (to the left) was found, then we are the
		// next biggest one.
		if (retVal == NULL) {
		    return self;
		} else {
		    return retVal;
		}
    } else {
		// The current element is strictly smaller than the one specified.
		// We have to find a bigger one.
		if (self->getRight() != NULL) {
		    return findFirstGreaterThanOrEqualTo(self->getRight(), key);
		} else {
		    return NULL;
		}
    }
}

template<class elementClass, class keyClass, class nodeClass> unsigned
AVLNodeLogic<elementClass,keyClass,nodeClass>::checkAndReturnDepth(
	nodeClass			self,
    unsigned			*countedElements)
{
    // We've been inserted and not deleted
    _ASSERT(self->getBalance() != AVLNew);

    (*countedElements)++;

    // Assert that the links all match up.
    _ASSERT(!self->getLeft() || self->getLeft()->getParent() == self);
    _ASSERT(!self->getRight() || self->getRight()->getParent() == self);

    // The basic binary tree ordering property applies
    _ASSERT(!self->getRight() || self->getKey() < self->getRight()->getKey());
    _ASSERT(!self->getLeft() || self->getKey() > self->getLeft()->getKey());

    // The AVL balance property applies
    unsigned leftDepth;
    if (self->getLeft()) {
		leftDepth = checkAndReturnDepth(self->getLeft(), countedElements);
    } else {
		leftDepth = 0;
    }

    unsigned rightDepth;
    if (self->getRight()) {
		rightDepth = checkAndReturnDepth(self->getRight(), countedElements);
    } else {
		rightDepth = 0;
    }

    if (leftDepth == rightDepth) {
        _ASSERT(self->getBalance() == AVLBalanced);
		return(leftDepth + 1);
    }

    if (leftDepth == rightDepth + 1) {
        _ASSERT(self->getBalance() == AVLLeft);
		return(leftDepth + 1);
    }

    if (leftDepth + 1 == rightDepth) {
        _ASSERT(self->getBalance() == AVLRight);
		return(rightDepth + 1);
    }

    _ASSERT(!"AVL Tree out of balance");
    return(0);
}

template<class elementClass, class keyClass, class nodeClass> void
AVLNodeLogic<elementClass,keyClass,nodeClass>::insert(
	nodeClass							self,
    AVLTree<elementClass,keyClass>		*intoTree,
    elementClass						*element)
{
    _ASSERT(intoTree);
    _ASSERT(self->getLeft() == NULL && self->getRight() == NULL && self->getParent() == NULL);

	self->setElement(element);
    _ASSERT(self->getElement());

    intoTree->insertions++;

    // Special case the empty tree case.
    if (intoTree->tree == NULL) {
		intoTree->tree = self;
		self->setBalance(AVLBalanced);
		// We already know all of the links are NULL, which is correct for this case.
		return;
    }

    // Find the leaf position at which to do this insertion.

    nodeClass currentNode = intoTree->tree;
	//
	// There is no need to set this, since it's never referenced before its set, but PREfast can't figure that out, 
	// so I'm changing the code to support the tools.  Yay!
	//
    nodeClass previousNode = NULL;
    while (currentNode) {
		previousNode = currentNode;
		if (currentNode->getElement()->getKey() < element->getKey()) {
		    currentNode = currentNode->getRight();
		} else {
            _ASSERT((currentNode->getElement()->getKey() > element->getKey()) && "Trying to insert a duplicate item.  Use something other than an AVL tree.");
		    currentNode = currentNode->getLeft();
		}
    }

    self->setBalance(AVLBalanced);
    self->setParent(previousNode);
    _ASSERT(self->getParent());

    if (previousNode->getElement()->getKey() < self->getElement()->getKey()) {
        _ASSERT(!previousNode->getRight());
		previousNode->setRight(self);
		rightAdded(previousNode, intoTree);
//		intoTree->check();
    } else { 
        _ASSERT(!previousNode->getLeft());
		previousNode->setLeft(self);
		leftAdded(previousNode, intoTree);
//		intoTree->check();
    }
}

template<class elementClass, class keyClass, class nodeClass> void
AVLNodeLogic<elementClass,keyClass,nodeClass>::rightAdded(
	nodeClass						self,
    AVLTree<elementClass,keyClass>	*tree)
{
    //We've just gotten one deeper on our right side.
    _ASSERT(self->getBalance() != AVLNew);
    
    if (self->getBalance() == AVLLeft) {
		self->setBalance(AVLBalanced);
		// The depth of the subtree rooted here hasn't changed, we're done
		return;
    }

    if (self->getBalance() == AVLBalanced) {
		// We've just gotten one deeper, but are still balanced.  Update and recurse up the
		// tree.
		self->setBalance(AVLRight);
		if (self->getParent()) {
		    if (self->getParent()->getRight() == self) {
				rightAdded(self->getParent(), tree);
		    } else {
                _ASSERT(self->getParent()->getLeft() == self);
				leftAdded(self->getParent(), tree);
		    }
		}
		return;
    }

    _ASSERT(self->getBalance() == AVLRight);
    // We've just gone to double right (ie, out of balance).
    _ASSERT(self->getRight());
    if (self->getRight()->getBalance() == AVLRight) {
		singleRotate(self,tree,self->getRight(),AVLRight);
    } else {
        _ASSERT(self->getRight()->getBalance() == AVLLeft);	// Else we shouldn't have been AVLRight before the call
		doubleRotate(self,tree,self->getRight(),self->getRight()->getLeft(),AVLRight);
    }
}

template<class elementClass, class keyClass, class nodeClass> void
AVLNodeLogic<elementClass,keyClass,nodeClass>::leftAdded(
	nodeClass						self,
    AVLTree<elementClass,keyClass>	*tree)
{
    //We've just gotten one deeper on our right side.
    _ASSERT(self->getBalance() != AVLNew);
    
    if (self->getBalance() == AVLRight) {
		self->setBalance(AVLBalanced);
		// The depth of the subtree rooted here hasn't changed, we're done
		return;
    }

    if (self->getBalance() == AVLBalanced) {
		// We've just gotten one deeper, but are still balanced.  Update and recurse up the
		// tree.
		self->setBalance(AVLLeft);
		if (self->getParent()) {
		    if (self->getParent()->getRight() == self) {
				rightAdded(self->getParent(), tree);
		    } else {
                _ASSERT(self->getParent()->getLeft() == self);
				leftAdded(self->getParent(), tree);
		    }
		}
		return;
    }

    _ASSERT(self->getBalance() == AVLLeft);
    // We've just gone to double left (ie, out of balance).
    _ASSERT(self->getLeft());
    if (self->getLeft()->getBalance() == AVLLeft) {
		singleRotate(self,tree,self->getLeft(),AVLLeft);
    } else {
        _ASSERT(self->getLeft()->getBalance() == AVLRight);	// Else we shouldn't have been AVLLeft before the call
		doubleRotate(self,tree,self->getLeft(),self->getLeft()->getRight(),AVLLeft);
    }
}

template<class elementClass, class keyClass, class nodeClass> void
AVLNodeLogic<elementClass,keyClass,nodeClass>::singleRotate(
	nodeClass								self,
    AVLTree<elementClass,keyClass>			*tree,
    nodeClass								child,
    AVLBalance								 whichSide)
{
    // self is the parent node.

    _ASSERT(tree);
    _ASSERT(child);
    _ASSERT(whichSide == AVLRight || whichSide == AVLLeft);

    _ASSERT(whichSide != AVLRight || self->getRight() == child);
    _ASSERT(whichSide != AVLLeft || self->getLeft() == child);

    tree->singleRotations++;

    // Promote the child to our position in the tree.

    if (self->getParent()) {
		if (self->getParent()->getLeft() == self) {
		    self->getParent()->setLeft(child);
		    child->setParent(self->getParent());
		} else {
            _ASSERT(self->getParent()->getRight() == self);
		    self->getParent()->setRight(child);
		    child->setParent(self->getParent());
		}
    } else {
		// We're the root of the tree
        _ASSERT(tree->tree == self);
		tree->tree = child;
		child->setParent(NULL);
    }

    // Attach the child's light subtree to our heavy side (ie., where the child is attached now)
    // Then, attach us to the child's light subtree
    if (whichSide == AVLRight) {
		self->setRight(child->getLeft());
		if (self->getRight()) {
		    self->getRight()->setParent(self);
		}

		child->setLeft(self);
		self->setParent(child);
    } else {
		self->setLeft(child->getRight());
		if (self->getLeft()) {
		    self->getLeft()->setParent(self);
		}

		child->setRight(self);
		self->setParent(child);	
    }

    // Finally, now both our and our (former) child's balance is "balanced"
    self->setBalance(AVLBalanced);
    child->setBalance(AVLBalanced);
    // NB. One of the cases in delete will result in the above balance settings being incorrect.  That
    // case fixes up the settings after we return.
}

template<class elementClass, class keyClass, class nodeClass> void
AVLNodeLogic<elementClass,keyClass,nodeClass>::doubleRotate(
	nodeClass								self,
    AVLTree<elementClass,keyClass>			*tree,
    nodeClass								child,
    nodeClass								grandchild,
    AVLBalance								 whichSide)
{
    _ASSERT(tree && child && grandchild);
    _ASSERT(whichSide == AVLLeft || whichSide == AVLRight);

    _ASSERT(whichSide != AVLLeft || (self->getLeft() == child && child->getBalance() == AVLRight));
    _ASSERT(whichSide != AVLRight || (self->getRight() == child && child->getBalance() == AVLLeft));

    _ASSERT(child->getParent() == self);
    _ASSERT(grandchild->getParent() == child);

    tree->doubleRotations++;

    // Write down a copy of all of the subtrees; see Knuth v3 p454 for the picture.
    // NOTE: The alpha and delta trees are never moved, so we don't store them.
    nodeClass beta;
    nodeClass gamma;

    if (whichSide == AVLRight) {
		beta = grandchild->getLeft();
		gamma = grandchild->getRight();
    } else {
		beta = grandchild->getRight();
		gamma = grandchild->getLeft();
    }

    // Promote grandchild to our position
    if (self->getParent()) {
        if (self->getParent()->getLeft() == self) {
		    self->getParent()->setLeft(grandchild);
        } else {
            _ASSERT(self->getParent()->getRight() == self);
		    self->getParent()->setRight(grandchild);
    	}
    } else {
        _ASSERT(tree->tree == self);
		tree->tree = grandchild;
    }
    grandchild->setParent(self->getParent());

    // Attach the appropriate children to grandchild
    if (whichSide == AVLRight) {
		grandchild->setRight(child);
		grandchild->setLeft(self);
    } else {
		grandchild->setRight(self);
		grandchild->setLeft(child);
    }
    self->setParent(grandchild);
    child->setParent(grandchild);

    // Attach beta and gamma to us and child.
    if (whichSide == AVLRight) {
		self->setRight(beta);
		if (beta) {
		    beta->setParent(self);
		}
		child->setLeft(gamma);
		if (gamma) {
		    gamma->setParent(child);
		}
    } else {
		self->setLeft(beta);
		if (beta) {
	    	beta->setParent(self);
		}
		child->setRight(gamma);
		if (gamma) {
		    gamma->setParent(child);
		}
    }

    // Now update the balance fields.
    switch (grandchild->getBalance()) {
	case AVLLeft:
		if (whichSide == AVLRight) {
		    self->setBalance(AVLBalanced);
		    child->setBalance(AVLRight);
		} else {
		    self->setBalance(AVLRight);
		    child->setBalance(AVLBalanced);
		}
		break;

	case  AVLBalanced:
		self->setBalance(AVLBalanced);
		child->setBalance(AVLBalanced);
		break;

	case AVLRight:
		if (whichSide == AVLRight) {
		    self->setBalance(AVLLeft);
		    child->setBalance(AVLBalanced);
		} else {
		    self->setBalance(AVLBalanced);
		    child->setBalance(AVLLeft);
		}
		break;

	default:
        _ASSERT(!"Bogus balance value");
    }
    grandchild->setBalance(AVLBalanced);
}

template<class elementClass, class keyClass, class nodeClass> void
AVLNodeLogic<elementClass,keyClass,nodeClass>::remove(
	nodeClass								self,
    AVLTree<elementClass,keyClass>			*fromTree)
{
    _ASSERT(fromTree);
    _ASSERT(self->getBalance() == AVLRight || self->getBalance() == AVLLeft || self->getBalance() == AVLBalanced);

    fromTree->deletions++;

    if (self->getLeft() == NULL) {
		// The right child either doesn't exist or is a leaf (because of the AVL balance property)
        _ASSERT((!self->getRight() && self->getBalance() == AVLBalanced) ||
	       (self->getBalance() == AVLRight
			&& self->getRight()->getBalance() == AVLBalanced
			&& self->getRight()->getRight() == NULL
			&& self->getRight()->getLeft() == NULL));

		if (self->getRight()) {
		    self->getRight()->setParent(self->getParent());
		}

		if (self->getParent()) {
		    if (self->getParent()->getLeft() == self) {
				self->getParent()->setLeft(self->getRight());
		        gotOneShorter(self->getParent(),fromTree,AVLLeft);
		    } else {
                _ASSERT(self->getParent()->getRight() == self);
				self->getParent()->setRight(self->getRight());
				gotOneShorter(self->getParent(),fromTree,AVLRight);
		    }
		} else {
            _ASSERT(fromTree->tree == self);
		    fromTree->tree = self->getRight();
		}
    } else if (self->getRight() == NULL) {
		// The left child must be a left because of the AVL balance property
        _ASSERT(self->getLeft()
			&& self->getBalance() == AVLLeft
			&& self->getLeft()->getBalance() == AVLBalanced
			&& self->getLeft()->getRight() == NULL
			&& self->getLeft()->getLeft() == NULL);
		self->getLeft()->setParent(self->getParent());
		if (self->getParent()) {
		    if (self->getParent()->getLeft() == self) {
				self->getParent()->setLeft(self->getLeft());
		        gotOneShorter(self->getParent(),fromTree,AVLLeft);
		    } else {
                _ASSERT(self->getParent()->getRight() == self);
				self->getParent()->setRight(self->getLeft());
				gotOneShorter(self->getParent(),fromTree,AVLRight);
		    }
		} else {
            _ASSERT(fromTree->tree == self);
		    fromTree->tree = self->getLeft();
		}
    } else {
		// Find the symmetric successor and promote it.  The symmetric successor is the smallest element in the right
		// subtree; it's found by following all left links in the right subtree until we find a node with no left link.
		// That node may be promoted to the place of this without corrupting the binary tree ordering properties. (We could
		// just as easily use the symmetric predecessor by finding the largest element in the left subtree, but there's
		// no point.)

		nodeClass successorCandidate = self->getRight();
		while (successorCandidate->getLeft()) {
		    successorCandidate = successorCandidate->getLeft();
		}

		nodeClass shorterRoot;
		AVLBalance shorterSide;
		if (successorCandidate->getParent()->getLeft() == successorCandidate) {
		    // We need to promote the successor's child (if any) to its position, then
	    	// promote it to our position.
		    shorterRoot = successorCandidate->getParent();
		    shorterSide = AVLLeft;
	    	successorCandidate->getParent()->setLeft(successorCandidate->getRight());
		    if (successorCandidate->getRight()) {
				successorCandidate->getRight()->setParent(successorCandidate->getParent());
		    }
	    
		    successorCandidate->setRight(self->getRight());
		    successorCandidate->setLeft(self->getLeft());
	    	successorCandidate->setBalance(self->getBalance());
		    successorCandidate->getRight()->setParent(successorCandidate);
		    successorCandidate->getLeft()->setParent(successorCandidate);
	    	if (self->getParent()) {
				if (self->getParent()->getLeft() == self) {
				    self->getParent()->setLeft(successorCandidate);
				} else {
                    _ASSERT(self->getParent()->getRight() == self);
				    self->getParent()->setRight(successorCandidate);
				}
		    } else {
                _ASSERT(fromTree->tree == self);
				fromTree->tree = successorCandidate;
		    }
		    successorCandidate->setParent(self->getParent());
		} else {
		    // The successor was our child, just directly promote it.
            _ASSERT(successorCandidate->getParent() == self);
		    if (self->getParent()) {
	    	    if (self->getParent()->getRight() == self) {
				    self->getParent()->setRight(successorCandidate);
		        } else {
                    _ASSERT(self->getParent()->getLeft() == self);
				    self->getParent()->setLeft(successorCandidate);
		    	}
		    } else {
                _ASSERT(fromTree->tree == self);
				fromTree->tree = successorCandidate;
		    }
		    successorCandidate->setParent(self->getParent());
		    successorCandidate->setLeft(self->getLeft());
		    if (self->getLeft()) {
				self->getLeft()->setParent(successorCandidate);
		    }
		    // We just made our right subtree shorter.
		    successorCandidate->setBalance(self->getBalance());
		    shorterRoot = successorCandidate;
		    shorterSide = AVLRight;
		}
		if (shorterRoot) {
		    gotOneShorter(shorterRoot,fromTree,shorterSide);
		}
    }

    self->setBalance(AVLNew);
	self->setLeft(NULL);
	self->setRight(NULL);
	self->setParent(NULL);
	self->setElement(NULL);
//    fromTree->check();
}

template<class elementClass, class keyClass, class nodeClass> void
AVLNodeLogic<elementClass,keyClass,nodeClass>::gotOneShorter(
	nodeClass							self,
    AVLTree<elementClass,keyClass>		*tree,
    AVLBalance			 				whichSide)
{
    _ASSERT(whichSide == AVLLeft || whichSide == AVLRight);

    if (self->getBalance() == AVLBalanced) {
		// We've just shrunk one subttree, but our depth has stayed the same.
		// Reset our balance indicator and punt.
		if (whichSide == AVLRight) {
		    self->setBalance(AVLLeft);
		} else {
		    self->setBalance(AVLRight);
		}
		return;
    } else if (self->getBalance() == whichSide) {
		// We just shrunk our heavy side; set our balance to neutral and recurse up the tree
		self->setBalance(AVLBalanced);
		if (self->getParent()) {
		    if (self->getParent()->getRight() == self) {
				gotOneShorter(self->getParent(),tree,AVLRight);
		    } else {
                _ASSERT(self->getParent()->getLeft() == self);
				gotOneShorter(self->getParent(),tree,AVLLeft);
		    }
		} // else we were the root; we're done
		return;
    } else {
		// We've just gone out of balance.  Figure out a rotation to do.  This is almost like having added a
		// node to the opposide side, except that the opposite side might be balanced.
		AVLBalance heavySide;
		nodeClass heavyChild;
		nodeClass replacement;
		if (whichSide == AVLRight) {
		    heavySide = AVLLeft;
		    heavyChild = self->getLeft();
		} else {
		    heavySide = AVLRight;
	    	heavyChild = self->getRight();
		}
        _ASSERT(heavyChild);
		if (heavyChild->getBalance() == heavySide) {
		    // Typical single rotation case
		    singleRotate(self,tree,heavyChild,heavySide);
		    replacement = heavyChild;
		} else if (heavyChild->getBalance() == whichSide) {
		    // Typical double rotation case
	    	nodeClass grandchild;
		    if (heavySide == AVLRight) {
				grandchild = heavyChild->getLeft();
	    	} else {
				grandchild = heavyChild->getRight();
		    }
		    doubleRotate(self,tree,heavyChild,grandchild,heavySide);
		    replacement = grandchild;
		} else {
            _ASSERT(heavyChild->getBalance() == AVLBalanced);
		    singleRotate(self,tree,heavyChild,heavySide);
		    // singleRotate has incorrectly set the balances; reset them
	    	self->setBalance(heavySide);
		    heavyChild->setBalance(whichSide);
		    // Overall depth hasn't changed; we're done.
	    	return;
		}

        // NB: we have now changed position in the tree, so parent, right & left have changed!
		if (!replacement->getParent()) {
		    // We just promoted our replacement to the root; we be done
		    return;
		}

		if (replacement->getParent()->getRight() == replacement) {
		    gotOneShorter(replacement->getParent(), tree,AVLRight);
		} else {
            _ASSERT(replacement->getParent()->getLeft() == replacement);
	    	gotOneShorter(replacement->getParent(), tree,AVLLeft);
		}
    }
}

template<class elementClass, class keyClass, class nodeClass> int
AVLNodeLogic<elementClass,keyClass,nodeClass>::inTree(
	nodeClass							self)
{
    return(self->getBalance() != AVLNew);
}

