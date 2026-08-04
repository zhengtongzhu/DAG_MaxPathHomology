import multiprocessing as mp


DENSITIES = [
    round((index + 1) * 0.1, 1)
    for index in range(9)
]


class TimeExperiment:
    """Measure runtime only."""

    def __init__(
        self,
        Rep_num: int,
        nodes_per_layer: list[int],
        seed: int | None = None,
    ) -> None:
        self.Rep_num = Rep_num
        self.nodes_per_layer = nodes_per_layer
        self.seed = seed

    def run(self) -> None:
        from .experiment_functions import compare_stratified

        compare_stratified(
            Rep_num=self.Rep_num,
            nodes_per_layer=self.nodes_per_layer,
            densities=DENSITIES,
            seed=self.seed,
        )


class MemoryExperiment:
    """Measure computation-specific peak additional RSS."""

    def __init__(
        self,
        Rep_num: int,
        nodes_per_layer: list[int],
        seed: int | None = None,
        sample_interval: float = 0.001,
        process_timeout: float | None = None,
    ) -> None:
        self.Rep_num = Rep_num
        self.nodes_per_layer = nodes_per_layer
        self.seed = seed
        self.sample_interval = sample_interval
        self.process_timeout = process_timeout

    def run(self) -> None:
        from .experiment_functions import compare_stratified_memory

        compare_stratified_memory(
            Rep_num=self.Rep_num,
            nodes_per_layer=self.nodes_per_layer,
            densities=DENSITIES,
            seed=self.seed,
            sample_interval=self.sample_interval,
            process_timeout=self.process_timeout,
        )


def run_time_experiments() -> None:
    """Run runtime experiments only."""
    experiment_table = [
        {
            "Rep_num": 200,
            "nodes_per_layer": [10, 10],
            "seed": 1,
        },
        {
            "Rep_num": 200,
            "nodes_per_layer": [10, 10, 10],
            "seed": 1,
        },
        {
            "Rep_num": 200,
            "nodes_per_layer": [4, 10, 10, 10],
            "seed": 1,
        },
    ]

    for parameters in experiment_table:
        experiment = TimeExperiment(**parameters)
        experiment.run()


def run_memory_experiments() -> None:
    """Run computation-specific peak-memory experiments only."""
    experiment_table = [
        {
            "Rep_num": 10,
            "nodes_per_layer": [10, 10],
            "seed": 1,
        },
        {
            "Rep_num": 10,
            "nodes_per_layer": [10, 10, 10],
            "seed": 1,
        },
        {
            "Rep_num": 10,
            "nodes_per_layer": [4, 10, 10, 10],
            "seed": 1,
        },
    ]

    for parameters in experiment_table:
        experiment = MemoryExperiment(
            **parameters,
            sample_interval=0.001,
            process_timeout=None,
        )
        experiment.run()

if __name__ == "__main__":
    mp.freeze_support()
    run_time = False
    run_memory = True
    
    if run_time:
        run_time_experiments()
    
    if run_memory:
        run_memory_experiments()