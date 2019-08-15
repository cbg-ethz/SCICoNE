//
// Created by Tuncel  Mustafa Anil on 2019-08-15.
// Based on user 992460's answer on stackoverflow.com/questions/8152720/
//

#ifndef SC_DNA_CUSTOMEXCEPTIONS_H
#define SC_DNA_CUSTOMEXCEPTIONS_H

#include <stdexcept>

#endif //SC_DNA_CUSTOMEXCEPTIONS_H

class NotImplementedException: public std::exception
{
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
    explicit NotImplementedException(const std::string& message):
            msg_(message)
    {}

    /** Destructor.
     * Virtual to allow for subclassing.
     */
    virtual ~NotImplementedException() throw (){}

    /** Returns a pointer to the (constant) error description.
     *  @return A pointer to a const char*. The underlying memory
     *          is in posession of the Exception object. Callers must
     *          not attempt to free the memory.
     */
    virtual const char* what() const throw (){
        return msg_.c_str();
    }

protected:
    /** Error message.
     */
    std::string msg_;
};