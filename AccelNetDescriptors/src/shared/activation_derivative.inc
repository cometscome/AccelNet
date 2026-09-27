        case(ACTIVATION_LINEAR); value = 1.0_real64
        case(ACTIVATION_TANH); value = 1.0_real64 - y*y
        case(ACTIVATION_LOGISTIC); value = y*(1.0_real64 - y)
        case(ACTIVATION_AENET_MTANH)
            tanhbx = tanh(b*x)
            value = a*b*(1.0_real64 - tanhbx*tanhbx)
        case(ACTIVATION_AENET_TWIST)
            tanhbx = tanh(b*x)
            value = a*b*(1.0_real64 - tanhbx*tanhbx) + c
        case(ACTIVATION_SOFTPLUS); value = 1.0_real64 - exp(-y)
        case(ACTIVATION_RELU); value = merge(1.0_real64, 0.0_real64, y > 0.0_real64)
        case(ACTIVATION_GAUSSIAN); value = -x*y
        case(ACTIVATION_COSINE); value = -sin(x)
        case(ACTIVATION_REVERSE_LOGISTIC); value = -y*(1.0_real64 - y)
        case(ACTIVATION_EXPONENTIAL); value = -y
        case(ACTIVATION_HARMONIC); value = 2.0_real64*x
