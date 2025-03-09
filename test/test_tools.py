import subprocess


progs = "../exe/"

create_matrix = progs + "CreateMatrix"
devide_probs = progs + "DevideUni"
prior_probs = progs + "PriorC"
merge_probs = progs + "MergePriorProbs"
em = progs + "Multi2Uni"


answ_err = "answ.err"
answ_out = "answ.out"
null_file = "null"

res_dir_corr = "corr_res/"
res_dir = "res/"
std_dir = "corr_std/"
inp = "input/"


def files_equal(name1: str, name2: str) -> bool:
	if name1.endswith("/") or name2.endswith("/"):
		return True
	f1 = open(name1, "r")
	f2 = open(name2, "r")

	assert f1.read() == f2.read()

	f1.close()
	f2.close()

	return True


def run(cmd, corr_err, corr_out, res, corr_res):
        with open(answ_out, 'w') as sout, open(answ_err, 'w') as serr:
                subprocess.Popen(cmd, stderr=serr, stdout=sout).communicate()

        files_equal(answ_err, std_dir+corr_err)
        files_equal(answ_out, std_dir+corr_out)
        files_equal(res_dir+res, res_dir_corr+corr_res)
