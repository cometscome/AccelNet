#include "CutoffFunction.h"
#include "NeuralNetwork.h"

extern "C" void n2p2_cutoff_value_derivative(int kind,
                                                double radius,
                                                double alpha,
                                                double distance,
                                                double* value,
                                                double* derivative)
{
    nnp::CutoffFunction cutoff;
    cutoff.setCutoffType(static_cast<nnp::CutoffFunction::CutoffType>(kind));
    cutoff.setCutoffRadius(radius);
    cutoff.setCutoffParameter(alpha);
    cutoff.fdf(distance, *value, *derivative);
}

extern "C" void n2p2_activation_value_derivative(int accelnetCode,
                                                   double input,
                                                   double* value,
                                                   double* derivative)
{
    using AF = nnp::NeuralNetwork::ActivationFunction;
    AF activation = AF::AF_UNSET;
    switch (accelnetCode)
    {
        case 0:  activation = AF::AF_IDENTITY;    break;
        case 1:  activation = AF::AF_TANH;        break;
        case 2:  activation = AF::AF_LOGISTIC;    break;
        case 5:  activation = AF::AF_RELU;        break;
        case 6:  activation = AF::AF_GAUSSIAN;    break;
        case 7:  activation = AF::AF_COS;         break;
        case 8:  activation = AF::AF_REVLOGISTIC; break;
        case 9:  activation = AF::AF_EXP;         break;
        case 10: activation = AF::AF_HARMONIC;    break;
        case 11: activation = AF::AF_SOFTPLUS;    break;
    }

    int const nodes[3] = {1, 1, 1};
    AF const activations[3] = {AF::AF_IDENTITY, activation, AF::AF_IDENTITY};
    double const connections[4] = {1.0, 0.0, 1.0, 0.0};
    nnp::NeuralNetwork network(3, nodes, activations);
    network.setConnections(connections);
    network.setInput(0, input);
    network.propagate();
    network.getOutput(value);
    network.calculateDEdG(derivative);
}
