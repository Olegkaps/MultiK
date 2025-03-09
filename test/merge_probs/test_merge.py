import pytest
import sys

sys.path.append("../")
from test_tools import *


@pytest.mark.parametrize(
	"cmd, corr_err, corr_out, res, corr_res",
	[
		([merge_probs],
			"no_args.err", "help.out", "", null_file),
		([merge_probs, "-e"],
			null_file, "help.out", "", null_file),
		([merge_probs, "-1", "abc", "-2", "cab", "-o"],
			null_file, "no_out.out", "", null_file),
		([merge_probs, "-o", "abc", "-2", "cab", "-1"],
			null_file, "no_1.out", "", null_file),
		([merge_probs, "-1", "abc", "-o", "cab", "-2"],
			null_file, "no_2.out", "", null_file),
		([merge_probs, "-1", "input/spline_down.txt", "-2", "input/spline_up.txt", "-o", "nodir/nofile"],
			"bad_out.err", "bad_out.out", "", null_file),
		([merge_probs, "-1", "input/nofile", "-2", "input/spline_up.txt", "-o", "res/nofile"],
			"bad_12.err", "bad_1.out", "", null_file),
		([merge_probs, "-1", "input/spline_down.txt", "-2", "input/nofile", "-o", "res/nofile"],
			"bad_12.err", "bad_2.out", "", null_file),
		([merge_probs, "--first", "input/spline_down.txt", "--second", "input/spline_up.txt", "--output", "res/file_50_1"],
			null_file, "", "file_50_1", "spline_merged.txt")
		# already non relevant
		#([merge_probs, "-1", "input/spline_down_longer.txt", "-2", "input/spline_up.txt", "-o", "res/firstbad"],
		#	"long_first.err", "long_first.out", "", null_file),
		#([merge_probs, "-2", "input/spline_down_longer.txt", "-1", "input/spline_up.txt", "-o", "res/secondbad"],
		#	"long_second.err", "long_second.out", "", null_file)

	]
) #spline_down_cap1.txt  spline_down_cap2.txt  - fully rewrite
def test_all_merge(cmd, corr_err, corr_out, res, corr_res):
	run(cmd, corr_err, corr_out, res, corr_res)
