using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;


// Unit test for the ASELib AVLTree class.

namespace AVLTest
{
    class Program
    {

        static void CheckBoolResult(bool expectedResult, bool actualResult, string errorMessage)
        {
            if (expectedResult != actualResult)
            {
                throw new ArgumentException(errorMessage);
            }
        }

        class TestTreeElement : IComparable
        {
            public TestTreeElement(int key_)
            {
                key = key_;
                value = key * 100000;
            }

            public int key;
            public int value = 0;
            public bool inTree = false;

            public int CompareTo(object genericPeer)
            {
                var peer = (TestTreeElement)genericPeer;
                return key.CompareTo(peer.key);
            }
        }

        static void Main(string[] args)
        {
            var fortyTwo = new TestTreeElement(42);
            var random = new Random(12345); // Use a specified seed in order to get deterministic behavior at least at first.
            for (int iteration = 0; iteration < 100; iteration++)
            {
                var tree = new ASETools.AVLTree<TestTreeElement>();
                tree.check();

                var operationLog = new List<string>();

                TestTreeElement returnedValue;

                CheckBoolResult(false, tree.findMin(out returnedValue), "Tree should be empty");
                CheckBoolResult(false, tree.findMax(out returnedValue), "Tree should be empty");
                CheckBoolResult(false, tree.Lookup(fortyTwo, out returnedValue), "Tree should be empty");

                int arraySize = 10000;
                var array = new TestTreeElement[arraySize];

                for (int i = 0; i < arraySize; i++)
                {
                    array[i] = new TestTreeElement(i);
                }

                int nInTree = 0;

                for (int innerIteration = 0; innerIteration < arraySize * 3; innerIteration++)
                {
                    int whichElement = random.Next(arraySize);

                    CheckBoolResult(array[whichElement].inTree, tree.Lookup(array[whichElement], out returnedValue), "Tree should contain element iff array does");

                    if (array[whichElement].inTree)
                    {
                        CheckBoolResult(true, array[whichElement].value == returnedValue.value, "Returned value from tree is incorrect");
                        CheckBoolResult(true, tree[array[whichElement]].value == array[whichElement].value, "Indexing into tree should return the correct value");

                        if (random.Next(2) == 0)
                        {
                            //
                            // Overwrite it.
                            //
                            array[whichElement].value++;
                            operationLog.Add("Overwrite at " + whichElement + " with " + array[whichElement].value);
                            tree[array[whichElement]].value = array[whichElement].value;
                        } else
                        {
                            //
                            // Delete it.
                            //
                            operationLog.Add("Delete " + whichElement);
                            tree.Delete(array[whichElement]);
                            array[whichElement].inTree = false;
                            nInTree--;
                        }
                    } else
                    {
                        //
                        // Insert it into the tree, half the time directly an half the time using an assignment.
                        //
                        if (random.Next(2) == 0)
                        {
                            array[whichElement].value++;
                            operationLog.Add("Insert " + array[whichElement].value + " at " + whichElement);


                            tree.Insert(array[whichElement]);
                        }
                        else
                        {
                            operationLog.Add("Insert by assignment " + array[whichElement].value + " at " + whichElement);
                            tree[array[whichElement]] = array[whichElement];
                        }

                        array[whichElement].inTree = true;
                        nInTree++;

                        tree.check();
                    }
                } // inner iteration


                int greatestValueAt = -1;
                int smallestValueAt = -1;
                for (int arrayIndex = 0; arrayIndex < arraySize; arrayIndex++)
                {
                    if (array[arrayIndex].inTree)
                    {
                        if (smallestValueAt == -1)
                        {
                            smallestValueAt = arrayIndex;
                        }

                        greatestValueAt = arrayIndex;
                    }
                }

                CheckBoolResult(true, (nInTree != 0) == tree.findMin(out returnedValue) && (nInTree == 0 || returnedValue.key == smallestValueAt), "FindMin didn't work properly.");
                CheckBoolResult(true, (nInTree != 0) == tree.findMax(out returnedValue) && (nInTree == 0 || returnedValue.key == greatestValueAt), "FindMax didn't work properly.");

                //
                // Now walk the tree and find if it has the expected values.
                //
                for (int checkIteration = 0; checkIteration < arraySize; checkIteration++)
                {
                    CheckBoolResult(array[checkIteration].inTree, tree.Lookup(array[checkIteration], out returnedValue), "Element should be in tree iff it's marked that way in array.");
                    CheckBoolResult(true, !array[checkIteration].inTree || returnedValue.value == array[checkIteration].value, "Looked up incorrect element");
                    CheckBoolResult(true, !array[checkIteration].inTree || returnedValue.value == tree[array[checkIteration]].value, "Lokoed up incorrect element using indexing.");

                    //
                    // Should also check that FindFirstLessThanOrEqualTo returns the right result when the key isn't in the tree.
                    //

                    CheckBoolResult(true, !array[checkIteration].inTree || tree.FindFirstLessThanOrEqualTo(array[checkIteration], out returnedValue) && returnedValue.value == array[checkIteration].value, "FindFirstLessThanOrEqualTo didn't find the correct result");
                    CheckBoolResult(true, checkIteration <= smallestValueAt || smallestValueAt == -1 || tree.FindFirstLessThan(array[checkIteration], out returnedValue) && returnedValue.key < checkIteration, "FindFirstLessThan didn't work properly.");
                    CheckBoolResult(true, checkIteration <= smallestValueAt || array[returnedValue.key].inTree, "FindFirstLessThan returned an element not in the tree.");
                    if (checkIteration > smallestValueAt)
                    {
                        for (int arrayIndex = returnedValue.key + 1; arrayIndex < checkIteration; arrayIndex++)
                        {
                            CheckBoolResult(false, array[arrayIndex].inTree, "FindFirstLessThan seems to have missed something.");
                        }
                    }

                    CheckBoolResult(true, !array[checkIteration].inTree || tree.FindFirstGreaterThanOrEqualTo(array[checkIteration], out returnedValue) && returnedValue.value == array[checkIteration].value, "FindFirstGreaterThanOrEqualTo didn't find the correct result");
                    CheckBoolResult(true, checkIteration >= greatestValueAt || greatestValueAt == -1 || tree.FindFirstGreaterThan(array[checkIteration], out returnedValue) && returnedValue.key > checkIteration, "findFirstGreaterThan didn't work properly.");
                    CheckBoolResult(true, checkIteration >= greatestValueAt || array[returnedValue.key].inTree, "FindFirstLessThan returned an element not in the tree.");
                    if (checkIteration < greatestValueAt)
                    {
                        for (int arrayIndex = returnedValue.key + 1; arrayIndex < checkIteration; arrayIndex++)
                        {
                            CheckBoolResult(false, array[arrayIndex].inTree, "FindFirstLessThan seems to have missed something.");
                        }
                    }


                } // checkIteration


            } // iteration
        } // Main
    }
}
