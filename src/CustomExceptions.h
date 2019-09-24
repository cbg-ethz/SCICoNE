//
// Created by Tuncel  Mustafa Anil on 2019-08-15.
//

#ifndef SC_DNA_CUSTOMEXCEPTIONS_H
#define SC_DNA_CUSTOMEXCEPTIONS_H

#include <stdexcept>

#endif //SC_DNA_CUSTOMEXCEPTIONS_H

class InvalidMove: public std::exception
{
    /*
     * Gets thrown when the move is invalid given the state of the input tree
     * */
public:
    explicit InvalidMove(const char* message):
            msg_(message) {}

    explicit InvalidMove(std::string message):
            msg_(std::move(message)) {}

    ~InvalidMove() noexcept override=default;

    const char* what() const noexcept override{
        return msg_.c_str();
    }
protected:
    std::string msg_;
};

class InvalidTree: public std::exception
{
    /*
     * Gets thrown when an invalid tree is produced
     * */
public:
    explicit InvalidTree(const char* message):
            msg_(message) {}

    explicit InvalidTree(std::string message):
            msg_(std::move(message)) {}

    ~InvalidTree() noexcept override=default;

    const char* what() const noexcept override{
        return msg_.c_str();
    }
protected:
    std::string msg_;
};

class NotImplementedException: public std::exception
{
    // Based on user 992460's answer on stackoverflow.com/questions/8152720/

public:
    /** Constructor (C strings).
     *  @param message C-style string error message.
     *                 The string contents are copied upon construction.
     *                 Hence, responsibility for deleting the char* lies
     *                 with the caller.
     */
    explicit NotImplementedException(const char* message):
            msg_(message)
    {
    }

    /** Constructor (C++ STL strings).
     *  @param message The error message.
     */
    explicit NotImplementedException(std::string  message):
            msg_(std::move(message))
    {}

    /** Destructor.
     * Virtual to allow for subclassing.
     */
    ~NotImplementedException() noexcept override= default;

    /** Returns a pointer to the (constant) error description.
     *  @return A pointer to a const char*. The underlying memory
     *          is in posession of the Exception object. Callers must
     *          not attempt to free the memory.
     */
    const char* what() const noexcept override{
        return msg_.c_str();
    }

protected:
    /** Error message.
     */
    std::string msg_;
};

