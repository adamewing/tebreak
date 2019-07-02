import unittest

class MyTestCase(unittest.TestCase):

  def test_equal_output(self):
    expected = open("example.tab.reference.txt")
    output = open("example.tab.txt")
    
    self.assertEqual(len(expected.readlines()), len(output.readlines()))

    expected.close()
    output.close()
    
if __name__ == '__main__':
    unittest.main()
