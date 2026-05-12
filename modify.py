    def _run_param_mc(self, state):
        with open(f'Logs/alpha{self._device_id}.csv', 'a', newline='') as f:
            writer = csv.writer(f)
            writer.writerow([state.alpha])
        if not self._parameter_manager.has_parameters():
            return state

        if self._options.param_mcmc_steps is None:
            raise RuntimeError(
                "There are sampled parameters, but param_mcmc_steps is not set."
            )

        # energy = self.get_energy(state)
        # group_energy = self.get_group_energy(state)
        # l_pre = group_energy[-1]
        ri = 0
        step = 0
        ri1 = []
        l_pre1 = []
        n_res1 = []
        for _ in range(self._options.param_mcmc_steps):
            trial_params = self._parameter_manager.sample(state.parameters)
            if not self._parameter_manager.is_valid(trial_params):
                accept = False
            else:
                trial_state = SystemState(
                    state.positions,
                    state.velocities,
                    state.alpha,
                    state.energy,
                    state.group_energies,
                    state.box_vector,
                    trial_params,
                    state.mappings,
                )
                # trial_energy = self.get_energy(trial_state)
                trial_group = self.get_group_energies(trial_state)
                group = self.get_group_energies(state)
                y1 = trial_params.discrete.data
                y0 = state.parameters.discrete.data
                my_res = trial_group[0] + trial_group[1] + trial_group[3]
                l_res = group[0] + group[1] + group[3]
                l1 = ri * trial_group[-1]
                l0 = ri * group[-1]
                mse1 = trial_group[1] / sum(y1)
                mse0 = group[1] / sum(y0)
                # print(trial_group)
                l_pre1.append(l1)
                n_res1.append(sum(y1))
                # print(group)
                # print(y1, y0)
                # print(y1-y0)
                # print((sum(y1)-sum(y0))*(trial_group[1]/sum(y1)-group[1]/sum(y0)))
                # ri
                # print(np.log10(mse0/mse1)/(sum(y1)-sum(y0)))

                # delta = trial_energy - energy
                delta = (my_res + l1) - (l_res + l0)

                if delta < 0:
                    accept = True
                else:
                    if random.random() < math.exp(-delta):
                        accept = True
                    else:
                        accept = False

            if accept:
                state = trial_state
                step = step + 1
                delta = sum(y1) - sum(y0)
                # energy = trial_energy
                k1 = min(sum(y1), sum(y0))
                if delta != 0:
                    # ri = min(sum(y0) * max(np.log(mse0/mse1)/delta, 0), 10)
                    # ri = min(10 * max(np.log(mse0/mse1)/delta, 0), 10)
                    ri = 1 + min(max(k1 * np.log(mse0 / mse1) / delta, -1), 10)
                    ri1.append(ri)
                    # print(ri)
        with open(f'Logs/ri{self._device_id}.csv', 'a', newline='') as f:
            writer = csv.writer(f)
            writer.writerows([ri1])
        with open(f'Logs/step{self._device_id}.csv', 'a', newline='') as f:
            writer = csv.writer(f)
            writer.writerow([step])
        with open(f'Logs/Nres{self._device_id}.csv', 'a', newline='') as f:
            writer = csv.writer(f)
            writer.writerows([n_res1])
        with open(f'Logs/lpre{self._device_id}.csv', 'a', newline='') as f:
            writer = csv.writer(f)
            writer.writerows([l_pre1])
        # Update transfomers in case we rejected the
        # last MCMC move
        if not accept:
            self._transformers_update(state)

        return state
