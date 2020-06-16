import random


def run_simulate_computing():
    tmp_in_ = "tmp_out"
    tmp_out_ = "tmp_inp"

    N = 0
    with open(tmp_in_, 'r') as inf:
        with open(tmp_out_, 'w') as outf:
            for _ in inf.readlines():
                v = random.random() * 1000
                print(f"{v}", file=outf)
                N += 1
    print(f"Computed {N} proteins")


if __name__ == "__main__":
    run_simulate_computing()
