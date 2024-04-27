import cProfile
import unittest
import test_astar
import pstats

suite = unittest.TestLoader().loadTestsFromModule(test_astar)
cProfile.run("unittest.TextTestRunner().run(suite)", "profile_output.txt")

stats = pstats.Stats("profile_output.txt")
stats.sort_stats("cumulative")
stats.print_stats()
