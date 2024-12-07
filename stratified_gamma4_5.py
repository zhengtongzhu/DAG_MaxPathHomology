from experiment_func import recursive_stratified

class Experiment2:
    def __init__(self, Rep_num: int, nodes_per_layer: list[int], seed: int = None) -> None:
        self.Rep_num = Rep_num
        self.nodes_per_layer = nodes_per_layer
        self.seed = seed

    def run(self) -> None:
        densities = [round((i + 1) * 0.1, 1) for i in range(9)]
        return recursive_stratified(self.Rep_num, self.nodes_per_layer, densities, self.seed)

if __name__ == '__main__':

    experiment_table = [
    {"Rep_num": 200, "nodes_per_layer": [4, 10, 10, 10, 10], "seed": 1},
    {"Rep_num": 200, "nodes_per_layer": [4, 10, 10, 10, 10, 5], "seed": 1}
    ]

    for param in experiment_table:
        experiment = Experiment2(**param)
        experiment.run()