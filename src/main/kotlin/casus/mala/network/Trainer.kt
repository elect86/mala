package casus.mala.network

import casus.mala.common.Parameters
import casus.mala.dataHandling.DataHandler

/** A class for training a neural network. */
class Trainer(parametersFull: Parameters,
              network: Network,
              data: DataHandler) : Runner(parametersFull, network, data) {

    var finalTestLoss = Float.POSITIVE_INFINITY
    var initialTestLoss = Float.POSITIVE_INFINITY
    var finalValidationLoss = Float.POSITIVE_INFINITY
    var initialValidationLoss = Float.POSITIVE_INFINITY
//    self.optimizer = None
//    self.scheduler = None
    var patienceCounter = 0
    var lastEpoch = 0
//    self.last_loss = None
//    self.training_data_loaders = []
//    self.validation_data_loaders = []
//    self.test_data_loaders = []

//    # Samplers for the horovod case.
//    self.train_sampler = None
//    self.test_sampler = None
//    self.validation_sampler = None
//
//    self.__prepare_to_train(optimizer_dict)
//
//    self.tensor_board = None
//    self.full_visualization_path = None
//    if self.parameters.visualisation:
//    if not os.path.exists(self.parameters.visualisation_dir):
//    os.makedirs(self.parameters.visualisation_dir)
//    if self.parameters.visualisation_dir_append_date:
//    date_time = datetime.now().strftime("%Y-%m-%d-%H-%M-%S")
//    self.full_visualization_path = os.path.join(
//    self.parameters.visualisation_dir, date_time
//    )
//    os.makedirs(self.full_visualization_path)
//    else:
//    self.full_visualization_path = (
//    self.parameters.visualisation_dir
//    )
//
//    # Set the path to log files
//    self.tensor_board = SummaryWriter(self.full_visualization_path)
//    printout(
//    "Writing visualization output to",
//    self.full_visualization_path,
//    min_verbosity=1,
//    )
//
//    self.gradscaler = None
//    if self.parameters.use_mixed_precision:
//    printout("Using mixed precision via AMP.", min_verbosity=1)
//    self.gradscaler = torch.cuda.amp.GradScaler()
//
//    self.train_graph = None
//    self.validation_graph = None

    /** Train a network using data given by a DataHandler. */
    fun trainNetwork() {
        //###########################
        // CALCULATE INITIAL METRICS
        //###########################

        var tloss = Float.POSITIVE_INFINITY
        val vloss = validateNetwork("validation", parameters.afterBeforeTrainingMetric)

//        if self.data.test_data_sets:
//        tloss = self.__validate_network(
//            self.network,
//            "test",
//            self.parameters.after_before_training_metric,
//        )
//
//        # Collect and average all the losses from all the devices
//            if self.parameters_full.use_horovod:
//        vloss = self.__average_validation(vloss, "average_loss")
//        self.initial_validation_loss = vloss
//        if self.data.test_data_set is not None:
//        tloss = self.__average_validation(tloss, "average_loss")
//        self.initial_test_loss = tloss
//
//        printout(
//            "Initial Guess - validation data loss: ", vloss, min_verbosity=1
//        )
//        if self.data.test_data_sets:
//        printout(
//            "Initial Guess - test data loss: ", tloss, min_verbosity=1
//        )
//
//        # Save losses for later use.
//        self.initial_validation_loss = vloss
//        self.initial_test_loss = tloss
//
//        # Initialize all the counters.
//        checkpoint_counter = 0
//
//        # If we restarted from a checkpoint, we have to differently initialize
//        # the loss.
//        if self.last_loss is None:
//        vloss_old = vloss
//        else:
//        vloss_old = self.last_loss
//
//        ############################
//        # PERFORM TRAINING
//        ############################
//
//        for epoch in range(self.last_epoch, self.parameters.max_number_epochs):
//        start_time = time.time()
//
//        # Prepare model for training.
//    self.network.train()
//
//        # Process each mini batch and save the training loss.
//    training_loss_sum = torch.zeros(
//        1, device=self.parameters._configuration["device"]
//    )
//
//        # train sampler
//            if self.parameters_full.use_horovod:
//        self.train_sampler.set_epoch(epoch)
//
//        # shuffle dataset if necessary
//        if isinstance(self.data.training_data_sets[0], FastTensorDataset):
//        self.data.training_data_sets[0].shuffle()
//
//        if self.parameters._configuration["gpu"]:
//        torch.cuda.synchronize(
//            self.parameters._configuration["device"]
//        )
//        tsample = time.time()
//        t0 = time.time()
//        batchid = 0
//        for loader in self.training_data_loaders:
//        for inputs, outputs in loader:
//
//        if self.parameters.profiler_range is not None:
//        if batchid == self.parameters.profiler_range[0]:
//        torch.cuda.profiler.start()
//        if batchid == self.parameters.profiler_range[1]:
//        torch.cuda.profiler.stop()
//
//        torch.cuda.nvtx.range_push(f"step {batchid}")
//
//        torch.cuda.nvtx.range_push("data copy in")
//        inputs = inputs.to(
//            self.parameters._configuration["device"],
//            non_blocking=True,
//        )
//        outputs = outputs.to(
//            self.parameters._configuration["device"],
//            non_blocking=True,
//        )
//        # data copy in
//        torch.cuda.nvtx.range_pop()
//
//        loss = self.__process_mini_batch(
//            self.network, inputs, outputs
//        )
//        # step
//        torch.cuda.nvtx.range_pop()
//        training_loss_sum += loss
//
//        if (
//            batchid != 0
//            and (batchid + 1)
//            % self.parameters.training_report_frequency
//            == 0
//        ):
//        torch.cuda.synchronize(
//            self.parameters._configuration["device"]
//        )
//        sample_time = time.time() - tsample
//        avg_sample_time = (
//                sample_time
//                / self.parameters.training_report_frequency
//                          )
//        avg_sample_tput = (
//                self.parameters.training_report_frequency
//                * inputs.shape[0]
//                / sample_time
//                          )
//        printout(
//            f"batch {batchid + 1}, "  # /{total_samples}, "
//        f"train avg time: {avg_sample_time} "
//        f"train avg throughput: {avg_sample_tput}",
//        min_verbosity=2,
//        )
//        tsample = time.time()
//        batchid += 1
//        torch.cuda.synchronize(
//            self.parameters._configuration["device"]
//        )
//        t1 = time.time()
//        printout(f"training time: {t1 - t0}", min_verbosity=2)
//
//        training_loss = training_loss_sum.item() / batchid
//
//        # Calculate the validation loss. and output it.
//        torch.cuda.synchronize(
//            self.parameters._configuration["device"]
//        )
//        else:
//        batchid = 0
//        for loader in self.training_data_loaders:
//        for inputs, outputs in loader:
//        inputs = inputs.to(
//            self.parameters._configuration["device"]
//        )
//        outputs = outputs.to(
//            self.parameters._configuration["device"]
//        )
//        training_loss_sum += self.__process_mini_batch(
//            self.network, inputs, outputs
//        )
//        batchid += 1
//        training_loss = training_loss_sum.item() / batchid
//
//        vloss = self.__validate_network(
//            self.network,
//            "validation",
//            self.parameters.during_training_metric,
//        )
//
//        if self.parameters_full.use_horovod:
//        vloss = self.__average_validation(vloss, "average_loss")
//        if self.parameters_full.verbosity > 1:
//        printout(
//            "Epoch {0}: validation data loss: {1}, "
//            "training data loss: {2}".format(
//                epoch, vloss, training_loss
//            ),
//            min_verbosity=2,
//        )
//        else:
//        printout(
//            "Epoch {0}: validation data loss: {1}".format(
//                epoch, vloss
//            ),
//            min_verbosity=1,
//        )
//
//        # summary_writer tensor board
//        if self.parameters.visualisation:
//        self.tensor_board.add_scalars(
//            "Loss",
//            {"validation": vloss, "training": training_loss},
//            epoch,
//        )
//        self.tensor_board.add_scalar(
//            "Learning rate", self.parameters.learning_rate, epoch
//        )
//        if self.parameters.visualisation == 2:
//        for name, param in self.network.named_parameters():
//        self.tensor_board.add_histogram(name, param, epoch)
//        self.tensor_board.add_histogram(
//            f"{name}.grad", param.grad, epoch
//        )
//
//        # method to make sure that all pending events have been written
//        # to disk
//            self.tensor_board.close()
//
//        if self.parameters._configuration["gpu"]:
//        torch.cuda.synchronize(
//            self.parameters._configuration["device"]
//        )
//
//        # Mix the DataSets up (this function only does something
//        # in the lazy loading case).
//        if self.parameters.use_shuffling_for_samplers:
//        self.data.mix_datasets()
//        if self.parameters._configuration["gpu"]:
//        torch.cuda.synchronize(
//            self.parameters._configuration["device"]
//        )
//
//        # If a scheduler is used, update it.
//        if self.scheduler is not None:
//        if (
//            self.parameters.learning_rate_scheduler
//            == "ReduceLROnPlateau"
//        ):
//        self.scheduler.step(vloss)
//
//        # If early stopping is used, check if we need to do something.
//        if self.parameters.early_stopping_epochs > 0:
//        if vloss < vloss_old * (
//                1.0 - self.parameters.early_stopping_threshold
//                               ):
//        self.patience_counter = 0
//        vloss_old = vloss
//        else:
//        self.patience_counter += 1
//        printout(
//            "Validation accuracy has not improved enough.",
//            min_verbosity=1,
//        )
//        if (
//            self.patience_counter
//            >= self.parameters.early_stopping_epochs
//        ):
//        printout(
//            "Stopping the training, validation "
//            "accuracy has not improved for",
//            self.patience_counter,
//            "epochs.",
//            min_verbosity=1,
//        )
//        self.last_epoch = epoch
//        break
//
//        # If checkpointing is enabled, we need to checkpoint.
//        if self.parameters.checkpoints_each_epoch != 0:
//        checkpoint_counter += 1
//        if (
//            checkpoint_counter
//            >= self.parameters.checkpoints_each_epoch
//        ):
//        printout("Checkpointing training.", min_verbosity=0)
//        self.last_epoch = epoch
//        self.last_loss = vloss_old
//        self.__create_training_checkpoint()
//        checkpoint_counter = 0
//
//        printout(
//            "Time for epoch[s]:",
//            time.time() - start_time,
//            min_verbosity=2,
//        )
//
//        ############################
//        # CALCULATE FINAL METRICS
//        ############################
//
//        if (
//            self.parameters.after_before_training_metric
//            != self.parameters.during_training_metric
//        ):
//        vloss = self.__validate_network(
//            self.network,
//            "validation",
//            self.parameters.after_before_training_metric,
//        )
//        if self.parameters_full.use_horovod:
//        vloss = self.__average_validation(vloss, "average_loss")
//
//        # Calculate final loss.
//    self.final_validation_loss = vloss
//        printout("Final validation data loss: ", vloss, min_verbosity=0)
//
//        tloss = float("inf")
//        if len(self.data.test_data_sets) > 0:
//        tloss = self.__validate_network(
//            self.network,
//            "test",
//            self.parameters.after_before_training_metric,
//        )
//        if self.parameters_full.use_horovod:
//        tloss = self.__average_validation(tloss, "average_loss")
//        printout("Final test data loss: ", tloss, min_verbosity=0)
//        self.final_test_loss = tloss
//
//        # Clean-up for pre-fetching lazy loading.
//    if self.data.parameters.use_lazy_loading_prefetch:
//        self.training_data_loaders.cleanup()
//        self.validation_data_loaders.cleanup()
//        if len(self.data.test_data_sets) > 0:
//        self.test_data_loaders.cleanup()
    }

    /** Validate a network, using test or validation data. */
    fun validateNetwork(dataSetType: String, validation: Parameters.TrainingMetric) {
        if (dataSetType == "test") {
            TODO()
//            data_loaders = self.test_data_loaders
//            data_sets = self.data.test_data_sets
//            number_of_snapshots = self.data.nr_test_snapshots
//            offset_snapshots = (
//                    self.data.nr_validation_snapshots
//                    + self.data.nr_training_snapshots
//                               )
        }
        else if (dataSetType == "validation") {
//            data_loaders = self.validation_data_loaders
//            data_sets = self.data.validation_data_sets
//            number_of_snapshots = self.data.nr_validation_snapshots
//            offset_snapshots = self.data.nr_training_snapshots
        }
        else
            error("Please select test or validation when using this function.")
//        network.eval()
//        if validation_type == "ldos":
//        validation_loss_sum = torch.zeros(
//            1, device=self.parameters._configuration["device"]
//        )
//        with torch.no_grad():
//        if self.parameters._configuration["gpu"]:
//        report_freq = self.parameters.training_report_frequency
//        torch.cuda.synchronize(
//            self.parameters._configuration["device"]
//        )
//        tsample = time.time()
//        batchid = 0
//        for loader in data_loaders:
//        for x, y in loader:
//        x = x.to(
//            self.parameters._configuration["device"],
//            non_blocking=True,
//        )
//        y = y.to(
//            self.parameters._configuration["device"],
//            non_blocking=True,
//        )
//
//        if (
//            self.parameters.use_graphs
//                    and self.validation_graph is None
//        ):
//        printout(
//            "Capturing CUDA graph for validation.",
//            min_verbosity=2,
//        )
//        s = torch.cuda.Stream(
//            self.parameters._configuration["device"]
//        )
//        s.wait_stream(
//            torch.cuda.current_stream(
//                self.parameters._configuration[
//                    "device"
//                ]
//            )
//        )
//        # Warmup for graphs
//        with torch.cuda.stream(s):
//        for _ in range(20):
//        with torch.cuda.amp.autocast(
//                enabled=self.parameters.use_mixed_precision
//                                    ):
//        prediction = network(x)
//        loss = network.calculate_loss(
//            prediction, y
//        )
//        torch.cuda.current_stream(
//            self.parameters._configuration["device"]
//        ).wait_stream(s)
//
//        # Create static entry point tensors to graph
//        self.static_input_validation = (
//                torch.empty_like(x)
//                                       )
//        self.static_target_validation = (
//                torch.empty_like(y)
//                                        )
//
//        # Capture graph
//            self.validation_graph = torch.cuda.CUDAGraph()
//        with torch.cuda.graph(self.validation_graph):
//        with torch.cuda.amp.autocast(
//                enabled=self.parameters.use_mixed_precision
//                                    ):
//        self.static_prediction_validation = (
//                network(
//                    self.static_input_validation
//                )
//                                            )
//        self.static_loss_validation = network.calculate_loss(
//            self.static_prediction_validation,
//            self.static_target_validation,
//        )
//
//        if self.validation_graph:
//        self.static_input_validation.copy_(x)
//        self.static_target_validation.copy_(y)
//        self.validation_graph.replay()
//        validation_loss_sum += (
//                self.static_loss_validation
//                               )
//        else:
//        with torch.cuda.amp.autocast(
//                enabled=self.parameters.use_mixed_precision
//                                    ):
//        prediction = network(x)
//        loss = network.calculate_loss(
//            prediction, y
//        )
//        validation_loss_sum += loss
//        if (
//            batchid != 0
//            and (batchid + 1) % report_freq == 0
//        ):
//        torch.cuda.synchronize(
//            self.parameters._configuration["device"]
//        )
//        sample_time = time.time() - tsample
//        avg_sample_time = sample_time / report_freq
//        avg_sample_tput = (
//                report_freq * x.shape[0] / sample_time
//                          )
//        printout(
//            f"batch {batchid + 1}, "  # /{total_samples}, "
//        f"validation avg time: {avg_sample_time} "
//        f"validation avg throughput: {avg_sample_tput}",
//        min_verbosity=2,
//        )
//        tsample = time.time()
//        batchid += 1
//        torch.cuda.synchronize(
//            self.parameters._configuration["device"]
//        )
//        else:
//        batchid = 0
//        for loader in data_loaders:
//        for x, y in loader:
//        x = x.to(self.parameters._configuration["device"])
//        y = y.to(self.parameters._configuration["device"])
//        prediction = network(x)
//        validation_loss_sum += network.calculate_loss(
//            prediction, y
//        ).item()
//        batchid += 1
//
//        validation_loss = validation_loss_sum.item() / batchid
//        return validation_loss
//        elif (
//            validation_type == "band_energy"
//                    or validation_type == "total_energy"
//        ):
//        errors = []
//        if isinstance(
//            self.validation_data_loaders, MultiLazyLoadDataLoader
//        ):
//        loader_id = 0
//        for loader in data_loaders:
//        grid_size = self.data.parameters.snapshot_directories_list[
//            loader_id + offset_snapshots
//        ].grid_size
//
//        actual_outputs = np.zeros(
//            (grid_size, self.data.output_dimension)
//        )
//        predicted_outputs = np.zeros(
//            (grid_size, self.data.output_dimension)
//        )
//        last_start = 0
//
//        for x, y in loader:
//
//        x = x.to(self.parameters._configuration["device"])
//        length = int(x.size()[0])
//        predicted_outputs[
//            last_start : last_start + length, :
//        ] = self.data.output_data_scaler.inverse_transform(
//        self.network(x).to("cpu"), as_numpy=True
//        )
//        actual_outputs[last_start : last_start + length, :] = (
//        self.data.output_data_scaler.inverse_transform(
//            y, as_numpy=True
//        )
//        )
//
//        last_start += length
//        errors.append(
//            self._calculate_energy_errors(
//                actual_outputs,
//                predicted_outputs,
//                validation_type,
//                loader_id + offset_snapshots,
//            )
//        )
//        loader_id += 1
//
//        else:
//        for snapshot_number in range(
//            offset_snapshots, number_of_snapshots + offset_snapshots
//        ):
//        # Get optimal batch size and number of batches per snapshotss
//            grid_size = self.data.parameters.snapshot_directories_list[
//        snapshot_number
//    ].grid_size
//
//        optimal_batch_size = self._correct_batch_size_for_testing(
//            grid_size, self.parameters.mini_batch_size
//        )
//        number_of_batches_per_snapshot = int(
//            grid_size / optimal_batch_size
//        )
//
//        actual_outputs, predicted_outputs = (
//        self._forward_entire_snapshot(
//            snapshot_number,
//            data_sets[0],
//            data_set_type[0:2],
//        number_of_batches_per_snapshot,
//        optimal_batch_size,
//        )
//        )
//
//        errors.append(
//            self._calculate_energy_errors(
//                actual_outputs,
//                predicted_outputs,
//                validation_type,
//                snapshot_number,
//            )
//        )
//        return np.mean(errors)
//        else:
//        raise Exception("Selected validation method not supported.")
    }
}