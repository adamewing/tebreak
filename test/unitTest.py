import unittest
import os

class MyTestCase(unittest.TestCase):

  def test_equal_output(self):
    expected = open("example.tab.reference.txt" , "rw+")
#    tail -n +2 example.tab.txt > example.tab2.txt 
#    sed 's/^.\{36\}//' example.tab2.txt > example.tab.txt
    output = open("example.tab1.txt" , "rw+")
    
    self.assertEquals(expected.readlines(), output.readlines())
    
if __name__ == '__main__':
    unittest.main()
