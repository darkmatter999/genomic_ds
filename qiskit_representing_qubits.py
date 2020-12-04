'''
Whereas in classical computing the bit is the 'basic variable', in quantum computing this is the 'qubit'. Both bit and qubit can only store
a single piece of information, they can only ever output a binary information, either 0 or 1, on or off.
However, qubits are different from classical bits in that they can be in a kind of 'in-between-state' (between 0 and 1) instead of just
being wholely in 0 or wholely in 1. They can be in whichever linear combination of states 0 and 1. In quantum mechanics this is called a 
'superposition'.

While classical bits have the binary value, 0 or 1, assigned as an integer to the variable, qubits store 'statevectors', 
orthogonal column vectors where [1,0] represents state 0 and [0,1] represents 1.

The bra-ket notation from quantum mechanics is used to represent these statevectors. A bra is a row vector and is initialized as <x|.
A ket represents a column vector and is initialized as |x>.

As long as the magnitude (np.linalg.norm) of the statevector is 1, i.e. the total probability of measuring the vector states sums to 1, states other
than 0 or 1 can be added to the statevector, such as:

|q0> = [1/sqrt(2), i/sqrt(2)]

There is usually a 'complex' (i, or j in Python) part in a statevector.

Code measuring the probability of the state being 1:

'''
# Import all required tools
from qiskit import QuantumCircuit, execute, Aer
from qiskit.visualization import plot_histogram, plot_bloch_vector
from math import sqrt, pi

qc = QuantumCircuit(1)  # Create a quantum circuit with one qubit
initial_state = [0,1]   # Define initial_state as |1>. The probability of this state will subsequently be measured.
qc.initialize(initial_state, 0) # Apply initialisation operation to the 0th qubit
qc.draw()  # Let's view our circuit

backend = Aer.get_backend('statevector_simulator') # A simulator (not real quantum machine) for statevectors is used.

result = execute(qc,backend).result() # Do the simulation, returning the result
out_state = result.get_statevector()
print(out_state) # Display the output state vector

'''
Here, as expected, the state vector [0,1] - however, with a complex part attached - is output.
'''

qc.measure_all() # Measuring in Qiskit means assigning a probability to what the output of a given qubit will be.
qc.draw()
counts = result.get_counts()
plot_histogram(counts)

'''
The resulting figure shows, as expected, a 100% probability of the output being |1>.

Now, instead of the state vector being [0,1] we use the 'quantum statevector' [1/sqrt(2), i/sqrt(2)]:
'''

initial_state = [1/sqrt(2), 1j/sqrt(2)]  # Define state |q_0>

qc = QuantumCircuit(1) # Must redefine qc. As in the previous example, the Quantum Circuit still consists of a single qubit.
qc.initialize(initial_state, 0) # Initialise the 0th qubit in the state `initial_state`
state = execute(qc,backend).result().get_statevector() # Execute the circuit
print(state)           # Print the result

results = execute(qc,backend).result().get_counts()
plot_histogram(results)

'''
The resulting figure shows that the probability of measuring |0> or |1> is equal, i.e. 50% for 0 and 50% for 1.
Why is this?
If we want to find the probability of measuring a state |x> in the state |psi> we use the following equation:
p(|x>) = |<x|psi>|² --> squared magnitude of the inner product of |x> and |psi>

Following above 'complex' statevector, we calculate step-by-step as follows:

|q0> = 1/sqrt(2)|0> + i/sqrt(2)|1> --> this is just our qubit state

<0|q0> = 1/sqrt(2)<0|0> + i/sqrt(2)<0|1> --> The probability of measuring the output 0 given above state |q0>.
                                             The <0|0> and <0|1> represent the inner products of |x> and |psi>.

       = 1/sqrt(2) * 1 + i/sqrt(2) * 0 --> The inner (dot) products of [0, 0], i.e. the statevectors [1,0] [1,0], and
                                           [0, 1], i.e. the statevectors [1,0] [0,1], are 1 and 0, respectively

       = 1/sqrt(2)

therefore --> |<0|q0>|² (squared magnitude of the inner product, following above formula p(|x>) = |<x|psi>|²) = 1/2

If we verified |1> instead of |0> we'd gotten the same result, since, as stated above, the probabilities are 50/50.

This above reasoning defines how information is gotten out of quantum states.

