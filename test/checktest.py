import unittest
import os

class MyTestCase(unittest.TestCase):

  def test_equal_output(self):
    expected = open("example.tab.reference.txt" , "rw+")
    output = open("example.tab.txt" , "rw+")
    
    self.assertEquals(len(expected.readlines()), len(output.readlines()))
    
if __name__ == '__main__':
    unittest.main()
