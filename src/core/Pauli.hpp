#pragma once

#include <string>
#include <vector>
#include <map>
#include <complex>
#include <algorithm>

namespace pmrqmc::core {

// Enum class for Pauli operators
enum class PauliOp : char {
    I = 'I',
    X = 'X', 
    Y = 'Y',
    Z = 'Z'
};

// Convert PauliOp to string
inline std::string to_string(PauliOp op) {
    return std::string(1, static_cast<char>(op));
}

// Represents a Pauli string as a product of Pauli operators on specific qubits
class PauliString {
private:
    std::map<int, PauliOp> operators_;
    
public:
    PauliString() = default;
    
    // Construct from a string like "X0 Z1 Y2"
    explicit PauliString(const std::string& pauli_str) {
        std::string cleaned = pauli_str;
        cleaned.erase(std::remove_if(cleaned.begin(), cleaned.end(), ::isspace), cleaned.end());
        
        for (size_t i = 0; i < cleaned.size(); ) {
            if (i + 1 >= cleaned.size()) break;
            
            char op_char = cleaned[i];
            int qubit = cleaned[i + 1] - '0';
            
            if (op_char == 'X') operators_[qubit] = PauliOp::X;
            else if (op_char == 'Y') operators_[qubit] = PauliOp::Y;
            else if (op_char == 'Z') operators_[qubit] = PauliOp::Z;
            // 'I' operators are implicit (not stored)
            
            i += 2;
        }
    }
    
    // Add an operator to a specific qubit
    void add_operator(int qubit, PauliOp op) {
        if (op != PauliOp::I) {
            operators_[qubit] = op;
        } else {
            // Remove I operator (identity is implicit)
            operators_.erase(qubit);
        }
    }
    
    // Remove operator from a specific qubit
    void remove_operator(int qubit) {
        operators_.erase(qubit);
    }
    
    // Get operator at specific qubit (returns I if not present)
    PauliOp get_operator(int qubit) const {
        auto it = operators_.find(qubit);
        if (it != operators_.end()) {
            return it->second;
        }
        return PauliOp::I;
    }
    
    // Check if operator exists at qubit (non-identity)
    bool has_operator(int qubit) const {
        return operators_.find(qubit) != operators_.end();
    }
    
    // Get all qubits with non-identity operators
    std::vector<int> get_qubits() const {
        std::vector<int> qubits;
        qubits.reserve(operators_.size());
        for (const auto& [qubit, op] : operators_) {
            qubits.push_back(qubit);
        }
        return qubits;
    }
    
    // Get number of non-identity operators
    size_t size() const {
        return operators_.size();
    }
    
    // Check if this is the identity operator
    bool is_identity() const {
        return operators_.empty();
    }
    
    // Multiply two Pauli strings (returns new PauliString)
    PauliString operator*(const PauliString& other) const {
        PauliString result = *this;
        
        for (const auto& [qubit, op] : other.operators_) {
            auto it = result.operators_.find(qubit);
            if (it == result.operators_.end()) {
                // No operator at this qubit, just add the new one
                if (op != PauliOp::I) {
                    result.operators_[qubit] = op;
                }
            } else {
                // Operator exists, multiply them
                PauliOp existing_op = it->second;
                PauliOp new_op = multiply_operators(existing_op, op);
                
                if (new_op == PauliOp::I) {
                    result.operators_.erase(it);
                } else {
                    it->second = new_op;
                }
            }
        }
        
        return result;
    }
    
    // In-place multiplication
    PauliString& operator*=(const PauliString& other) {
        *this = *this * other;
        return *this;
    }
    
    // Check if two Pauli strings commute
    bool commutes_with(const PauliString& other) const {
        // Two Pauli strings commute if they share an even number of anti-commuting pairs
        int anti_commuting_pairs = 0;
        
        // Get all qubits involved in either operator
        std::vector<int> all_qubits = this->get_qubits();
        std::vector<int> other_qubits = other.get_qubits();
        
        // Combine and remove duplicates
        all_qubits.insert(all_qubits.end(), other_qubits.begin(), other_qubits.end());
        std::sort(all_qubits.begin(), all_qubits.end());
        all_qubits.erase(std::unique(all_qubits.begin(), all_qubits.end()), all_qubits.end());
        
        for (int qubit : all_qubits) {
            PauliOp op1 = this->get_operator(qubit);
            PauliOp op2 = other.get_operator(qubit);
            
            if (operators_anti_commute(op1, op2)) {
                anti_commuting_pairs++;
            }
        }
        
        return anti_commuting_pairs % 2 == 0;
    }
    
    // Convert to string representation
    std::string to_string() const {
        if (operators_.empty()) {
            return "I";
        }
        
        std::string result;
        for (const auto& [qubit, op] : operators_) {
            result += static_cast<char>(op);
            result += std::to_string(qubit) + " ";
        }
        
        if (!result.empty()) {
            result.pop_back(); // Remove trailing space
        }
        
        return result;
    }
    
    // Apply to a basis state (bit string)
    // Returns the resulting state and phase factor
    std::pair<std::string, std::complex<double>> apply(const std::string& state) const {
        std::string result_state = state;
        std::complex<double> phase = 1.0;
        
        for (const auto& [qubit, op] : operators_) {
            if (qubit >= 0 && qubit < static_cast<int>(state.size())) {
                char bit = state[qubit];
                
                switch (op) {
                    case PauliOp::X: // Bit flip
                        result_state[qubit] = (bit == '0') ? '1' : '0';
                        break;
                    case PauliOp::Z: // Phase flip
                        if (bit == '1') {
                            phase *= -1.0;
                        }
                        break;
                    case PauliOp::Y: // Both bit flip and phase flip
                        result_state[qubit] = (bit == '0') ? '1' : '0';
                        if (bit == '1') {
                            phase *= std::complex<double>(0, -1.0); // -i
                        } else {
                            phase *= std::complex<double>(0, 1.0); // +i
                        }
                        break;
                    case PauliOp::I:
                        // Identity, do nothing
                        break;
                }
            }
        }
        
        return {result_state, phase};
    }
    
    // Equality comparison - order independent
    bool operator==(const PauliString& other) const {
        if (operators_.size() != other.operators_.size()) {
            return false;
        }
        
        for (const auto& [qubit, op] : operators_) {
            auto it = other.operators_.find(qubit);
            if (it == other.operators_.end() || it->second != op) {
                return false;
            }
        }
        
        return true;
    }
    
    bool operator!=(const PauliString& other) const {
        return !(*this == other);
    }
    
    // Get hash for use in unordered containers
    size_t hash() const {
        size_t seed = 0;
        for (const auto& [qubit, op] : operators_) {
            seed ^= std::hash<int>{}(qubit) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            seed ^= static_cast<size_t>(op) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }

private:
    // Multiply two Pauli operators
    static PauliOp multiply_operators(PauliOp op1, PauliOp op2) {
        // Pauli multiplication table (X*X = I, Y*Y = I, Z*Z = I, X*Y = Z, Y*Z = X, Z*X = Y, etc.)
        // Note: phases are ignored here, only operator part matters for our purposes
        if (op1 == PauliOp::I) return op2;
        if (op2 == PauliOp::I) return op1;
        if (op1 == op2) return PauliOp::I;  // X*X = I, Y*Y = I, Z*Z = I
        
        // Different operators
        if ((op1 == PauliOp::X && op2 == PauliOp::Y) || 
            (op1 == PauliOp::Y && op2 == PauliOp::X)) {
            return PauliOp::Z;
        }
        if ((op1 == PauliOp::Y && op2 == PauliOp::Z) || 
            (op1 == PauliOp::Z && op2 == PauliOp::Y)) {
            return PauliOp::X;
        }
        if ((op1 == PauliOp::Z && op2 == PauliOp::X) || 
            (op1 == PauliOp::X && op2 == PauliOp::Z)) {
            return PauliOp::Y;
        }
        
        return PauliOp::I;  // Shouldn't reach here
    }
    
    // Check if two operators anti-commute
    static bool operators_anti_commute(PauliOp op1, PauliOp op2) {
        if (op1 == PauliOp::I || op2 == PauliOp::I) return false;
        if (op1 == op2) return false;
        return true; // Any two different non-identity Pauli operators anti-commute
    }
};

// Hash function for PauliString
struct PauliStringHash {
    size_t operator()(const PauliString& p) const {
        return p.hash();
    }
};

}