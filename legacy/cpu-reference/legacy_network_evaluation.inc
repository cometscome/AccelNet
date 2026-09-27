! Standalone NN/reference evaluator retained from gpu checkpoint 543b180.
    subroutine evaluate_network(self, input, output)
        class(atomic_network), intent(in) :: self
        real(real64), intent(in) :: input(:)
        real(real64), intent(out) :: output
        real(real64) :: values(self%maxnodes, self%nlayers)
        integer :: layer, i, j, nin, nout, offset
        values = 0.0_real64
        values(1:self%nodes(1), 1) = input
        do layer = 1, self%nlayers - 1
            nin = self%nodes(layer); nout = self%nodes(layer + 1)
            offset = self%weight_offsets(layer) + 1
            do j = 1, nout
                values(j, layer + 1) = self%weights(offset + nin*nout + j - 1)
                do i = 1, nin
                    values(j, layer + 1) = values(j, layer + 1) + &
                        self%weights(offset + (i - 1)*nout + j - 1)*values(i, layer)
                end do
                values(j, layer + 1) = activate(values(j, layer + 1), self%activation(layer))
            end do
        end do
        output = values(1, self%nlayers)
    end subroutine evaluate_network

    subroutine network_input_gradient(self, input, output, gradient)
        class(atomic_network), intent(in) :: self
        real(real64), intent(in) :: input(:)
        real(real64), intent(out) :: output, gradient(:)
        real(real64) :: values(self%maxnodes, self%nlayers), preactivation(self%maxnodes, self%nlayers)
        real(real64) :: delta(self%maxnodes), previous(self%maxnodes)
        integer :: layer, i, j, nin, nout, offset
        values = 0.0_real64; preactivation = 0.0_real64; values(1:self%nodes(1), 1) = input
        do layer = 1, self%nlayers - 1
            nin = self%nodes(layer); nout = self%nodes(layer + 1); offset = self%weight_offsets(layer) + 1
            do j = 1, nout
                values(j, layer + 1) = self%weights(offset + nin*nout + j - 1)
                do i = 1, nin
                    values(j, layer + 1) = values(j, layer + 1) + self%weights(offset + (i - 1)*nout + j - 1)*values(i, layer)
                end do
                preactivation(j, layer + 1) = values(j, layer + 1)
                values(j, layer + 1) = activate(preactivation(j, layer + 1), self%activation(layer))
            end do
        end do
        output = values(1, self%nlayers); delta = 0.0_real64; delta(1) = 1.0_real64
        do layer = self%nlayers - 1, 1, -1
            nin = self%nodes(layer); nout = self%nodes(layer + 1); offset = self%weight_offsets(layer) + 1
            do j = 1, nout
                delta(j) = delta(j)*activation_derivative(preactivation(j, layer + 1), &
                    values(j, layer + 1), self%activation(layer))
            end do
            previous = 0.0_real64
            do i = 1, nin
                do j = 1, nout
                    previous(i) = previous(i) + self%weights(offset + (i - 1)*nout + j - 1)*delta(j)
                end do
            end do
            delta = previous
        end do
        gradient = delta(1:self%nodes(1))
    end subroutine network_input_gradient

