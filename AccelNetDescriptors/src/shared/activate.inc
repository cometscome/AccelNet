        case(ACTIVATION_LINEAR); y = x
        case(ACTIVATION_TANH); y = tanh(x)
        case(ACTIVATION_LOGISTIC); y = 1.0_real64/(1.0_real64 + exp(-x))
        case(ACTIVATION_AENET_MTANH)
            y = a*tanh(b*x)
        case(ACTIVATION_AENET_TWIST)
            y = a*tanh(b*x) + c*x
        case(ACTIVATION_SOFTPLUS)
            if (x > 0.0_real64) then
                y = x + log(1.0_real64 + exp(-x))
            else
                y = log(1.0_real64 + exp(x))
            end if
        case(ACTIVATION_RELU); y = max(x, 0.0_real64)
        case(ACTIVATION_GAUSSIAN); y = exp(-0.5_real64*x*x)
        case(ACTIVATION_COSINE); y = cos(x)
        case(ACTIVATION_REVERSE_LOGISTIC); y = 1.0_real64 - 1.0_real64/(1.0_real64 + exp(-x))
        case(ACTIVATION_EXPONENTIAL); y = exp(-x)
        case(ACTIVATION_HARMONIC); y = x*x
