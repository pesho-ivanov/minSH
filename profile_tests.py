import cProfile
import unittest
import astar_test
import pstats

suite = unittest.TestLoader().loadTestsFromModule(astar_test)
cProfile.run("unittest.TextTestRunner().run(suite)", "profile_output.txt")

stats = pstats.Stats("profile_output.txt")
stats.sort_stats("cumulative")
stats.print_stats()
