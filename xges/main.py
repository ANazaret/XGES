from evaluation.benchmarks import simulation
from xges.search import XGES


def main():
    data, dag = simulation(10, 1000, 1)

    model = XGES()
    model.fit(data)


if __name__ == '__main__':
    main()
