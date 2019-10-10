//
// CDDL HEADER START
//
// The contents of this file are subject to the terms of the Common Development
// and Distribution License Version 1.0 (the "License").
//
// You can obtain a copy of the license at
// http://www.opensource.org/licenses/CDDL-1.0.  See the License for the
// specific language governing permissions and limitations under the License.
//
// When distributing Covered Code, include this CDDL HEADER in each file and
// include the License file in a prominent location with the name LICENSE.CDDL.
// If applicable, add the following below this CDDL HEADER, with the fields
// enclosed by brackets "[]" replaced with your own identifying information:
//
// Portions Copyright (c) [yyyy] [name of copyright owner]. All rights reserved.
//
// CDDL HEADER END
//

//
// Copyright (c) 2019, Regents of the University of Minnesota.
// All rights reserved.
//
// Contributors:
//    Yaser Afshar
//


#ifndef SNAP_HELPER_HPP
#define SNAP_HELPER_HPP

/*!
 * \file helper.hpp
 *
 * \author Yaser Afshar (yafshar@umn.edu)
 *
 * \brief This file contains helper classes for memory management
 * which are STL like extension for multi-dimensional arrays
 *
 * \date 09-10-2019
 */

#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>

#ifdef HELPER_LOG_ERROR
#undef HELPER_LOG_ERROR
#endif

/*!
 * \brief Helper macro for printing error message
 * formats messages, filename, line number and function
 * name into an std::ostringstream object
 *
 */
#define HELPER_LOG_ERROR(msg)      \
  {                                \
    std::ostringstream ss;         \
    ss << "\nError :" << __FILE__  \
       << ":" << __LINE__ << ":@(" \
       << __FUNCTION__ << ")\n"    \
       << msg << "\n\n";           \
    std::cerr << ss.str();         \
  }

/*! Type alias for vector of constant dimension 3 */
using VectorOfSizeDIM = double[3];

/*! Type alias for vector of constant dimension 6 */
using VectorOfSizeSix = double[6];

/*! \class _Array_Basic The basic STL like container similar to `std::vector` to
 * handle multi-dimensional arrays.
 *
 * \brief An STL like container similar to <a
 * href="https://en.cppreference.com/w/cpp/container/vector">std::vector</a>
 * that encapsulates dynamic size arrays in a sequence container
 *
 * \tparam DataType The type of the elements. Default (double)
 */
template <class DataType>
class _Array_Basic
{
public:
  /*!
   * \brief Construct a new _Array_Basic object
   *
   */
  _Array_Basic();

  /*!
   * \brief Construct a new _Array_Basic object
   *
   * \param count The size of the container
   */
  explicit _Array_Basic(std::size_t const count);

  /*!
   * \brief Construct a new _Array_Basic object
   *
   * \param count The size of the container
   * \param value The value to initialize elements of the container with
   */
  _Array_Basic(std::size_t const count, DataType const value);

  /*!
   * \brief Construct a new _Array_Basic object
   *
   * \param count The size of the container
   * \param array The array of data to initialize elements of the container
   * with.
   */
  _Array_Basic(std::size_t const count, DataType const *array);

  /*!
   * \brief Construct a new _Array_Basic object
   * Copy constructor. Constructs the container with the copy of the contents of
   * other.
   *
   * \param other Another container to be used as source to initialize the
   * elements of the container with
   */
  _Array_Basic(_Array_Basic<DataType> const &other);

  /*!
   * \brief Construct a new _Array_Basic object
   * Move constructor. Constructs the container with the contents of other using
   * move semantics.
   *
   * \param other Another container to be used as source to initialize the
   * elements of the container with
   */
  _Array_Basic(_Array_Basic<DataType> &&other);

  /*!
   * \brief Destroy the _Array_Basic object
   *
   */
  ~_Array_Basic();

  /*!
   * \brief Copy assignment operator. Replaces the contents with a copy of the
   * contents of other
   *
   * \param other Another container to use as data source
   *
   * \return _Array_Basic<DataType>&
   */
  _Array_Basic<DataType> &operator=(_Array_Basic<DataType> const &other);

  /*!
   * \brief Move assignment operator. Replaces the contents with those of other
   * using move semantics
   *
   * \param other Another container to use as data source
   *
   * \return _Array_Basic<DataType>&
   */
  _Array_Basic<DataType> &operator=(_Array_Basic<DataType> &&other);

  /*!
   * \brief Returns pointer to the underlying array serving as element storage.
   *
   * \return DataType*
   */
  inline DataType const *data() const noexcept;
  inline DataType *data() noexcept;

  /*!
   * \brief Returns the number of elements in the container.
   *
   * \return std::size_t
   */
  inline std::size_t size() const;

  /*!
   * \brief Erases all elements from the container.
   *
   */
  inline void clear() noexcept;

  /*!
   * \brief Requests the removal of unused capacity.
   *
   */
  inline void shrink_to_fit();

  /*!
   * \brief Returns the number of elements that the container has currently
   * allocated space for.
   *
   * \return std::size_t
   */
  inline std::size_t capacity() const noexcept;

  /*!
   * \brief Appends the given element value to the end of the container.
   *
   * \param value
   */
  inline void push_back(DataType const &value);
  inline void push_back(DataType &&value);

protected:
  /*!
   * \brief Check the index range based on container size
   *
   * \param _n Index
   */
  inline void _range_check(int _n) const;

  /*!
   * \brief Check the index range based on input size
   *
   * \param _n Index
   * \param tsize Input size
   */
  inline void _range_check(int _n, std::size_t tsize) const;

protected:
  /*! Dynamic contiguous array */
  std::vector<DataType> m;
};

/*! \class Array1DView A 1-dimensional STL like container.
 *
 * \brief A 1-dimensional STL like container that encapsulates dynamic size
 * arrays in a sequence container
 *
 * \tparam DataType The type of the elements. Default (double)
 */
template <class DataType>
class Array1DView
{
public:
  /*!
   * \brief Construct a new Array1DView object
   *
   * \param count The size of the container
   * \param array The array of data to initialize elements of the container
   * with.
   */
  Array1DView(std::size_t const count, DataType *array);
  Array1DView(std::size_t const count, DataType const *array);

  /*!
   * \brief Construct a new Array1DView object
   * Copy constructor. Constructs the container with the copy of the contents of
   * other.
   *
   * \param other Another container to be used as source to initialize the
   * elements of the container with
   */
  Array1DView(Array1DView<DataType> const &other);

  /*!
   * \brief Destroy the Array1DView object
   *
   */
  ~Array1DView();

  /*!
   * \brief Returns pointer to the underlying array serving as element storage.
   *
   * \return DataType*
   */
  inline DataType const *data() const noexcept;
  inline DataType *data() noexcept;

  /*!
   * \brief Returns the element at specified location \b i.
   * No bounds checking is performed.
   *
   * \param i Position of the element to return
   *
   * \return const DataType The requested element.
   */
  inline const DataType operator()(int i) const;
  inline DataType &operator()(int i);

  /*!
   * \brief Returns the element at specified location \b i , with bounds
   * checking.
   *
   * \param i Position of the element to return
   *
   * \return const DataType The requested element.
   */
  inline DataType const at(int i) const;
  inline DataType &at(int i);

  /*!
   * \brief Returns the element at specified location \b i.
   * No bounds checking is performed.
   *
   * \param i Position of the element to return
   *
   * \return const DataType The requested element.
   */
  const DataType operator[](int i) const;
  DataType &operator[](int i);

private:
  Array1DView() = delete;

  Array1DView<DataType> &operator=(Array1DView<DataType> const &other) = delete;

  Array1DView<DataType> &operator=(Array1DView<DataType> &&other) = delete;

protected:
  /*!
   * \brief Check the index range based on input size
   *
   * \param _n Index
   * \param tsize Input size
   */
  inline void _range_check(int _n, std::size_t tsize) const;

protected:
  /*! The extent of the container in the 1st mode */
  std::size_t _extentZero;

  /*! Data pointer */
  DataType *const m;
};

template <class DataType>
class Array2DView
{
public:
  /*!
   * \brief Construct a new Array2DView object
   *
   * \param extentZero The extent of the container in the 1st mode
   * \param extentOne The extent of the container in the 2nd mode
   * \param array The array of data to set the pointer of the container to it.
   */
  Array2DView(std::size_t const extentZero,
              std::size_t const extentOne,
              DataType *array);

  Array2DView(std::size_t const extentZero,
              std::size_t const extentOne,
              DataType const *array);

  /*!
   * \brief Construct a new Array2DView object
   * Copy constructor. Constructs the container with the copy of the contents of
   * other.
   *
   * \param other Another container to be used as source to initialize the
   * elements of the container with
   */
  Array2DView(Array2DView<DataType> const &other);

  /*!
   * \brief Destroy the Array2D object
   *
   */
  ~Array2DView();

  /*!
   * \brief Returns pointer to the underlying array serving as element storage.
   *
   * \return DataType*
   */
  inline DataType const *data() const noexcept;
  inline DataType *data() noexcept;

  inline Array1DView<DataType> const data_1D(int i) const;
  inline Array1DView<DataType> data_1D(int i);

  /*!
   * \brief Returns the element at specified location \b (i, j).
   * No bounds checking is performed.
   *
   * \param i Position of the element in the 1st mode
   * \param j Position of the element in the 2nd mode
   *
   * \return const DataType The requested element.
   */
  inline const DataType operator()(int i, int j) const;
  inline DataType &operator()(int i, int j);

  /*!
   * \brief Returns the element at specified location \b (i, j) , with bounds
   * checking.
   *
   * \param i Position of the element in the 1st mode
   * \param j Position of the element in the 2nd mode
   *
   * \return const DataType The requested element.
   */
  inline DataType const at(int i, int j) const;
  inline DataType &at(int i, int j);

  /*! \class j_operator A helper class to provide multidimensional array access
   * semantics.
   *
   * \brief To provide 2-dimensional array access semantics, operator[] has to
   * return a reference to a 1D vector, which has to have its own operator[]
   * which returns a reference to the element.
   */
  class j_operator
  {
  public:
    /*!
     * \brief Construct a new j_operator object
     *
     * \param _array Refernce to Array2D class
     * \param i Position of the element in the 1st mode
     */
    j_operator(Array2DView<DataType> &_array, int i);

    /*!
     * \brief Provide array-like access and returns the element at specified
     * location \b [i][j]. No bounds checking is performed.
     *
     * \param j Position of the element in the 2nd mode
     *
     * \return const DataType The requested element.
     */
    const DataType operator[](int j) const;
    DataType &operator[](int j);

  private:
    /*! Refernce to Array2D class */
    Array2DView<DataType> &j_array;

    std::size_t _i;
  };

  /*!
   * \brief Provide array-like access and returns the element at specified
   * location \b [i][j]. No bounds checking is performed.
   *
   * \param i Position of the element in the 1st mode
   * \param j Position of the element in the 2nd mode
   *
   * \return const DataType The requested element.
   *
   * \note
   * To provide multidimensional array access semantics, we are using multiple
   * overloads for \code operator[] \endcode . For speed one should avoid this
   * complexity, uses \code operator() \endcode as \code (i, j) \endcode
   * directly.
   */
  const j_operator operator[](int i) const;
  j_operator operator[](int i);

private:
  Array2DView() = delete;

  Array2DView<DataType> &operator=(Array2DView<DataType> const &other) = delete;

  Array2DView<DataType> &operator=(Array2DView<DataType> &&other) = delete;

protected:
  /*!
   * \brief Check the index range based on input size
   *
   * \param _n Index
   * \param tsize Input size
   */
  inline void _range_check(int _n, std::size_t tsize) const;

protected:
  /*! The extent of the container in the 1st mode */
  std::size_t _extentZero;

  /*! The extent of the container in the 2nd mode */
  std::size_t _extentOne;

  /*! Data pointer */
  DataType *const m;
};

/*! \class Array2D A 2-dimensional STL like container.
 *
 * \brief A 2-dimensional STL like container that encapsulates dynamic size
 * arrays with a 2-dimensional shape in a sequence container
 *
 * \tparam DataType The type of the elements. Default (double)
 */
template <class DataType>
class Array2D : public _Array_Basic<DataType>
{
public:
  /*!
   * \brief Construct a new Array2D object
   *
   */
  Array2D();

  /*!
   * \brief Construct a new Array2D object
   *
   * \param extentZero The extent of the container in the 1st mode
   * \param extentOne The extent of the container in the 2nd mode
   */
  Array2D(std::size_t const extentZero, std::size_t const extentOne);

  /*!
   * \brief Construct a new Array2D object
   *
   * \param extentZero The extent of the container in the 1st mode
   * \param extentOne The extent of the container in the 2nd mode
   * \param value The value to initialize elements of the container with
   */
  Array2D(std::size_t const extentZero,
          std::size_t const extentOne,
          DataType const value);

  /*!
   * \brief Construct a new Array2D object
   *
   * \param extentZero The extent of the container in the 1st mode
   * \param extentOne The extent of the container in the 2nd mode
   * \param array The array of data to initialize elements of the container in a
   * row-major format.
   */
  Array2D(std::size_t const extentZero,
          std::size_t const extentOne,
          DataType const *array);

  /*!
   * \brief Construct a new Array2D object
   * Copy constructor. Constructs the container with the copy of the contents of
   * other.
   *
   * \param other Another container to be used as source to initialize the
   * elements of the container with
   */
  Array2D(Array2D<DataType> const &other);

  /*!
   * \brief Construct a new Array2D object
   * Move constructor. Constructs the container with the contents of other using
   * move semantics.
   *
   * \param other Another container to be used as source to initialize the
   * elements of the container with
   */
  Array2D(Array2D<DataType> &&other);

  /*!
   * \brief Destroy the Array2D object
   *
   */
  ~Array2D();

  /*!
   * \brief Copy assignment operator. Replaces the contents with a copy of the
   * contents of other
   *
   * \param other Another container to use as data source
   *
   * \return Array2D<DataType>&
   */
  Array2D<DataType> &operator=(Array2D<DataType> const &other);

  /*!
   * \brief Move assignment operator. Replaces the contents with those of other
   * using move semantics
   *
   * \param other Another container to use as data source
   *
   * \return Array2D<DataType>&
   */
  Array2D<DataType> &operator=(Array2D<DataType> &&other);

  /*!
   * \brief Returns Array1DView to the underlying starting element at row \c i.
   * \sa Array1DView
   *
   * \param i Row index
   * \return Array1DView<DataType>
   */
  inline Array1DView<DataType> const data_1D(int i) const;
  inline Array1DView<DataType> data_1D(int i);

  /*!
   * \brief Resizes the container to contain \c extentZero times \c extentOne
   * elements.
   *
   * \param extentZero
   * \param extentOne
   */
  inline void resize(int const extentZero, int const extentOne);

  /*!
   * \brief Resizes the container to contain \c extentZero times \c extentOne
   * elements.
   *
   * \param extentZero
   * \param extentOne
   * \param new_value The new value to initialize the new elements with
   */
  inline void resize(int const extentZero, int const extentOne, DataType const new_value);

  /*!
   * \brief Resizes the container to contain \c extentZero times \c extentOne
   * elements.
   *
   * \param extentZero
   * \param extentOne
   * \param new_array The new array of data to initialize elements of the
   * container with.
   */
  inline void resize(int const extentZero, int const extentOne, DataType const *new_array);

  /*!
   * \brief Returns the element at specified location \b (i, j).
   * No bounds checking is performed.
   *
   * \param i Position of the element in the 1st mode
   * \param j Position of the element in the 2nd mode
   *
   * \return const DataType The requested element.
   */
  inline const DataType operator()(int i, int j) const;
  inline DataType &operator()(int i, int j);

  /*!
   * \brief Returns the element at specified location \b (i, j) , with bounds
   * checking.
   *
   * \param i Position of the element in the 1st mode
   * \param j Position of the element in the 2nd mode
   *
   * \return const DataType The requested element.
   */
  inline DataType const at(int i, int j) const;
  inline DataType &at(int i, int j);

  /*! \class j_operator A helper class to provide multidimensional array access
   * semantics.
   *
   * \brief To provide 2-dimensional array access semantics, operator[] has to
   * return a reference to a 1D vector, which has to have its own operator[]
   * which returns a reference to the element.
   */
  class j_operator
  {
  public:
    /*!
     * \brief Construct a new j_operator object
     *
     * \param _array Refernce to Array2D class
     * \param i Position of the element in the 1st mode
     */
    j_operator(Array2D<DataType> &_array, int i);

    /*!
     * \brief Provide array-like access and returns the element at specified
     * location \b [i][j]. No bounds checking is performed.
     *
     * \param j Position of the element in the 2nd mode
     *
     * \return const DataType The requested element.
     */
    const DataType operator[](int j) const;
    DataType &operator[](int j);

  private:
    /*! Refernce to Array2D class */
    Array2D<DataType> &j_array;

    std::size_t _i;
  };

  /*!
   * \brief Provide array-like access and returns the element at specified
   * location \b [i][j]. No bounds checking is performed.
   *
   * \param i Position of the element in the 1st mode
   * \param j Position of the element in the 2nd mode
   *
   * \return const DataType The requested element.
   *
   * \note
   * To provide multidimensional array access semantics, we are using multiple
   * overloads for \code operator[] \endcode . For speed one should avoid this
   * complexity, uses \code operator() \endcode as \code (i, j) \endcode
   * directly.
   */
  const j_operator operator[](int i) const;
  j_operator operator[](int i);

protected:
  /*! The extent of the container in the 1st mode */
  std::size_t _extentZero;

  /*! The extent of the container in the 2nd mode */
  std::size_t _extentOne;
};

template <class DataType>
class Array3D : public _Array_Basic<DataType>
{
public:
  Array3D();

  Array3D(std::size_t const extentZero,
          std::size_t const extentOne,
          std::size_t const extentTwo);

  Array3D(std::size_t const extentZero,
          std::size_t const extentOne,
          std::size_t const extentTwo,
          DataType const value);

  Array3D(std::size_t const extentZero,
          std::size_t const extentOne,
          std::size_t const extentTwo,
          DataType const *array);

  Array3D(Array3D<DataType> const &other);

  Array3D(Array3D<DataType> &&other);

  ~Array3D();

  Array3D<DataType> &operator=(Array3D<DataType> const &other);

  Array3D<DataType> &operator=(Array3D<DataType> &&other);

  inline Array2DView<DataType> const data_2D(int i) const;
  inline Array2DView<DataType> data_2D(int i);

  inline Array1DView<DataType> const data_1D(int i, int j) const;
  inline Array1DView<DataType> data_1D(int i, int j);

  inline void resize(int const extentZero, int const extentOne, int const extentTwo);

  inline void resize(int const extentZero,
                     int const extentOne,
                     int const extentTwo,
                     DataType const new_value);

  inline void resize(int const extentZero,
                     int const extentOne,
                     int const extentTwo,
                     DataType const *new_array);

  inline const DataType operator()(int i, int j, int k) const;
  inline DataType &operator()(int i, int j, int k);

  inline DataType const at(int i, int j, int k) const;
  inline DataType &at(int i, int j, int k);

  class j_operator
  {
  public:
    j_operator(Array3D<DataType> &_array, int i);

    class k_operator
    {
    public:
      k_operator(Array3D<DataType> &_array, int i, int j);

      const DataType operator[](int k) const;
      DataType &operator[](int k);

    private:
      Array3D<DataType> &k_array;

      std::size_t _i;
      std::size_t _j;
    };

    const k_operator operator[](int j) const;
    k_operator operator[](int j);

  private:
    Array3D<DataType> &j_array;

    std::size_t _i;
  };

  const j_operator operator[](int i) const;
  j_operator operator[](int i);

protected:
  std::size_t _extentZero;
  std::size_t _extentOne;
  std::size_t _extentTwo;
};

template <class DataType>
_Array_Basic<DataType>::_Array_Basic() {}

template <class DataType>
_Array_Basic<DataType>::_Array_Basic(std::size_t const count) : m(count) {}

template <class DataType>
_Array_Basic<DataType>::_Array_Basic(std::size_t const count, DataType const value) : m(count, value) {}

template <class DataType>
_Array_Basic<DataType>::_Array_Basic(std::size_t const count, DataType const *array) : m(array, array + count) {}

template <class DataType>
_Array_Basic<DataType>::_Array_Basic(_Array_Basic<DataType> const &other) : m(other.m) {}

template <class DataType>
_Array_Basic<DataType>::_Array_Basic(_Array_Basic<DataType> &&other) : m(std::move(other.m)) {}

template <class DataType>
_Array_Basic<DataType>::~_Array_Basic() {}

template <class DataType>
_Array_Basic<DataType> &_Array_Basic<DataType>::operator=(_Array_Basic<DataType> const &other)
{
  m.resize(other.size());
  std::copy(other.m.begin(), other.m.end(), m.begin());
  return *this;
}

template <class DataType>
_Array_Basic<DataType> &_Array_Basic<DataType>::operator=(_Array_Basic<DataType> &&other)
{
  m = std::move(other.m);
  return *this;
}

template <class DataType>
inline DataType const *_Array_Basic<DataType>::data() const noexcept
{
  return m.data();
}

template <class DataType>
inline DataType *_Array_Basic<DataType>::data() noexcept
{
  return m.data();
}

template <class DataType>
inline std::size_t _Array_Basic<DataType>::size() const
{
  return m.size();
}

template <class DataType>
inline void _Array_Basic<DataType>::clear() noexcept
{
  m.clear();
}

template <class DataType>
inline void _Array_Basic<DataType>::shrink_to_fit()
{
  m.shrink_to_fit();
}

template <class DataType>
inline std::size_t _Array_Basic<DataType>::capacity() const noexcept
{
  return m.capacity();
}

template <class DataType>
inline void _Array_Basic<DataType>::push_back(DataType const &value)
{
  m.push_back(value);
}

template <class DataType>
inline void _Array_Basic<DataType>::push_back(DataType &&value)
{
  m.push_back(value);
}

template <class DataType>
inline void _Array_Basic<DataType>::_range_check(int _n) const
{
  if (_n >= size())
  {
    HELPER_LOG_ERROR("The input index is out of range! " + std::to_string(_n) + " >= " + std::to_string(size()));
    std::abort();
  }
}

template <class DataType>
inline void _Array_Basic<DataType>::_range_check(int _n, std::size_t tsize) const
{
  if (_n >= tsize)
  {
    HELPER_LOG_ERROR("The input index is out of range! " + std::to_string(_n) + " >= " + std::to_string(tsize));
    std::abort();
  }
}

template <class DataType>
Array1DView<DataType>::Array1DView(std::size_t const count, DataType *array) : _extentZero(count),
                                                                               m(array) {}

template <class DataType>
Array1DView<DataType>::Array1DView(std::size_t const count, DataType const *array) : _extentZero(count),
                                                                                     m(const_cast<DataType *>(array)) {}

template <class DataType>
Array1DView<DataType>::Array1DView(Array1DView<DataType> const &other) : _extentZero(other._extentZero),
                                                                         m(other.m) {}

template <class DataType>
Array1DView<DataType>::~Array1DView() {}

template <class DataType>
inline DataType const *Array1DView<DataType>::data() const noexcept
{
  return m;
}

template <class DataType>
inline DataType *Array1DView<DataType>::data() noexcept
{
  return m;
}

template <class DataType>
inline const DataType Array1DView<DataType>::operator()(int i) const
{
  return m[i];
}

template <class DataType>
inline DataType &Array1DView<DataType>::operator()(int i)
{
  return m[i];
}

template <class DataType>
inline DataType &Array1DView<DataType>::at(int i)
{
  _range_check(i, _extentZero);
  return m[i];
}

template <class DataType>
inline DataType const Array1DView<DataType>::at(int i) const
{
  _range_check(i, _extentZero);
  return m[i];
}

template <class DataType>
const DataType Array1DView<DataType>::operator[](int i) const
{
  return m[i];
}

template <class DataType>
DataType &Array1DView<DataType>::operator[](int i)
{
  return m[i];
}

template <class DataType>
inline void Array1DView<DataType>::_range_check(int _n, std::size_t tsize) const
{
  if (_n >= tsize)
  {
    HELPER_LOG_ERROR("The input index is out of range! " + std::to_string(_n) + " >= " + std::to_string(tsize));
    std::abort();
  }
}

template <class DataType>
Array2DView<DataType>::Array2DView(std::size_t const extentZero,
                                   std::size_t const extentOne,
                                   DataType *array) : _extentZero(extentZero),
                                                      _extentOne(extentOne),
                                                      m(array) {}

template <class DataType>
Array2DView<DataType>::Array2DView(std::size_t const extentZero,
                                   std::size_t const extentOne,
                                   DataType const *array) : _extentZero(extentZero),
                                                            _extentOne(extentOne),
                                                            m(const_cast<DataType *>(array)) {}

template <class DataType>
Array2DView<DataType>::Array2DView(Array2DView<DataType> const &other) : _extentZero(other._extentZero),
                                                                         _extentOne(other._extentOne),
                                                                         m(other.m) {}

template <class DataType>
Array2DView<DataType>::~Array2DView() {}

template <class DataType>
inline DataType const *Array2DView<DataType>::data() const noexcept
{
  return m;
}

template <class DataType>
inline DataType *Array2DView<DataType>::data() noexcept
{
  return m;
}

template <class DataType>
inline Array1DView<DataType> const Array2DView<DataType>::data_1D(int i) const
{
  return Array1DView<DataType>(_extentOne, m + i * _extentOne);
}

template <class DataType>
inline Array1DView<DataType> Array2DView<DataType>::data_1D(int i)
{
  return Array1DView<DataType>(_extentOne, m + i * _extentOne);
}

template <class DataType>
inline const DataType Array2DView<DataType>::operator()(int i, int j) const
{
  std::size_t const _n = i * _extentOne + j;
  return m[_n];
}

template <class DataType>
inline DataType &Array2DView<DataType>::operator()(int i, int j)
{
  std::size_t const _n = i * _extentOne + j;
  return m[_n];
}

template <class DataType>
inline DataType &Array2DView<DataType>::at(int i, int j)
{
  _range_check(i, _extentZero);
  _range_check(j, _extentOne);
  std::size_t const _n = i * _extentOne + j;
  return m[_n];
}

template <class DataType>
inline DataType const Array2DView<DataType>::at(int i, int j) const
{
  _range_check(i, _extentZero);
  _range_check(j, _extentOne);
  std::size_t const _n = i * _extentOne + j;
  return m[_n];
}

template <class DataType>
Array2DView<DataType>::j_operator::j_operator(Array2DView<DataType> &_array,
                                              int i) : j_array(_array),
                                                       _i(i) {}

template <class DataType>
const DataType Array2DView<DataType>::j_operator::operator[](int j) const
{
  std::size_t const _n = _i * j_array._extentOne + j;
  return j_array.m[_n];
}

template <class DataType>
DataType &Array2DView<DataType>::j_operator::operator[](int j)
{
  std::size_t const _n = _i * j_array._extentOne + j;
  return j_array.m[_n];
}

template <class DataType>
const typename Array2DView<DataType>::j_operator Array2DView<DataType>::operator[](int i) const
{
  return j_operator(*this, i);
}

template <class DataType>
typename Array2DView<DataType>::j_operator Array2DView<DataType>::operator[](int i)
{
  return j_operator(*this, i);
}

template <class DataType>
inline void Array2DView<DataType>::_range_check(int _n, std::size_t tsize) const
{
  if (_n >= tsize)
  {
    HELPER_LOG_ERROR("The input index is out of range! " + std::to_string(_n) + " >= " + std::to_string(tsize));
    std::abort();
  }
}

template <class DataType>
Array2D<DataType>::Array2D() : _Array_Basic<DataType>(),
                               _extentZero(0),
                               _extentOne(0) {}

template <class DataType>
Array2D<DataType>::Array2D(std::size_t const extentZero,
                           std::size_t const extentOne) : _Array_Basic<DataType>(extentZero * extentOne),
                                                          _extentZero(extentZero),
                                                          _extentOne(extentOne) {}

template <class DataType>
Array2D<DataType>::Array2D(std::size_t const extentZero,
                           std::size_t const extentOne,
                           DataType const value) : _Array_Basic<DataType>(extentZero * extentOne, value),
                                                   _extentZero(extentZero),
                                                   _extentOne(extentOne) {}

template <class DataType>
Array2D<DataType>::Array2D(std::size_t const extentZero,
                           std::size_t const extentOne,
                           DataType const *array) : _Array_Basic<DataType>(extentZero * extentOne, array),
                                                    _extentZero(extentZero),
                                                    _extentOne(extentOne) {}

template <class DataType>
Array2D<DataType>::Array2D(Array2D<DataType> const &other) : _Array_Basic<DataType>(other),
                                                             _extentZero(other._extentZero),
                                                             _extentOne(other._extentOne) {}

template <class DataType>
Array2D<DataType>::Array2D(Array2D<DataType> &&other) : _Array_Basic<DataType>(std::move(other)),
                                                        _extentZero(other._extentZero),
                                                        _extentOne(other._extentOne) {}

template <class DataType>
Array2D<DataType>::~Array2D() {}

template <class DataType>
Array2D<DataType> &Array2D<DataType>::operator=(Array2D<DataType> const &other)
{
  _Array_Basic<DataType>::operator=(other);
  _extentZero = other._extentZero;
  _extentOne = other._extentOne;
  return *this;
}

template <class DataType>
Array2D<DataType> &Array2D<DataType>::operator=(Array2D<DataType> &&other)
{
  _Array_Basic<DataType>::operator=(std::move(other));
  _extentZero = other._extentZero;
  _extentOne = other._extentOne;
  return *this;
}

template <class DataType>
inline Array1DView<DataType> const Array2D<DataType>::data_1D(int i) const
{
  return Array1DView<DataType>(_extentOne, this->m.data() + i * _extentOne);
}

template <class DataType>
inline Array1DView<DataType> Array2D<DataType>::data_1D(int i)
{
  return Array1DView<DataType>(_extentOne, this->m.data() + i * _extentOne);
}

template <class DataType>
inline void Array2D<DataType>::resize(int const extentZero, int const extentOne)
{
  _extentZero = extentZero;
  _extentOne = extentOne;
  std::size_t const _n = _extentZero * _extentOne;
  this->m.resize(_n);
}

template <class DataType>
inline void Array2D<DataType>::resize(int const extentZero,
                                      int const extentOne,
                                      DataType const new_value)
{
  _extentZero = extentZero;
  _extentOne = extentOne;
  std::size_t const _n = _extentZero * _extentOne;
  this->m.resize(_n, new_value);
}

template <class DataType>
inline void Array2D<DataType>::resize(int const extentZero,
                                      int const extentOne,
                                      DataType const *new_array)
{
  _extentZero = extentZero;
  _extentOne = extentOne;
  std::size_t const _n = _extentZero * _extentOne;
  this->m.resize(_n);
  std::copy(new_array, new_array + _n, this->m.data());
}

template <class DataType>
inline const DataType Array2D<DataType>::operator()(int i, int j) const
{
  std::size_t const _n = i * _extentOne + j;
  return this->m[_n];
}

template <class DataType>
inline DataType &Array2D<DataType>::operator()(int i, int j)
{
  std::size_t const _n = i * _extentOne + j;
  return this->m[_n];
}

template <class DataType>
inline DataType &Array2D<DataType>::at(int i, int j)
{
  this->_range_check(i, _extentZero);
  this->_range_check(j, _extentOne);
  std::size_t const _n = i * _extentOne + j;
  return this->m[_n];
}

template <class DataType>
inline DataType const Array2D<DataType>::at(int i, int j) const
{
  this->_range_check(i, _extentZero);
  this->_range_check(j, _extentOne);
  std::size_t const _n = i * _extentOne + j;
  return this->m[_n];
}

template <class DataType>
Array2D<DataType>::j_operator::j_operator(Array2D<DataType> &_array, int i) : j_array(_array),
                                                                              _i(i) {}

template <class DataType>
const DataType Array2D<DataType>::j_operator::operator[](int j) const
{
  std::size_t const _n = _i * j_array._extentOne + j;
  return j_array.m[_n];
}

template <class DataType>
DataType &Array2D<DataType>::j_operator::operator[](int j)
{
  std::size_t const _n = _i * j_array._extentOne + j;
  return j_array.m[_n];
}

template <class DataType>
const typename Array2D<DataType>::j_operator Array2D<DataType>::operator[](int i) const
{
  return j_operator(*this, i);
}

template <class DataType>
typename Array2D<DataType>::j_operator Array2D<DataType>::operator[](int i)
{
  return j_operator(*this, i);
}

template <class DataType>
Array3D<DataType>::Array3D() : _Array_Basic<DataType>(),
                               _extentZero(0),
                               _extentOne(0),
                               _extentTwo(0) {}

template <class DataType>
Array3D<DataType>::Array3D(std::size_t const extentZero,
                           std::size_t const extentOne,
                           std::size_t const extentTwo) : _Array_Basic<DataType>(extentZero * extentOne * extentTwo),
                                                          _extentZero(extentZero),
                                                          _extentOne(extentOne),
                                                          _extentTwo(extentTwo) {}

template <class DataType>
Array3D<DataType>::Array3D(std::size_t const extentZero,
                           std::size_t const extentOne,
                           std::size_t const extentTwo,
                           DataType const value) : _Array_Basic<DataType>(extentZero * extentOne * extentTwo, value),
                                                   _extentZero(extentZero),
                                                   _extentOne(extentOne),
                                                   _extentTwo(extentTwo) {}

template <class DataType>
Array3D<DataType>::Array3D(std::size_t const extentZero,
                           std::size_t const extentOne,
                           std::size_t const extentTwo,
                           DataType const *array) : _Array_Basic<DataType>(extentZero * extentOne * extentTwo, array),
                                                    _extentZero(extentZero),
                                                    _extentOne(extentOne),
                                                    _extentTwo(extentTwo) {}

template <class DataType>
Array3D<DataType>::Array3D(Array3D<DataType> const &other) : _Array_Basic<DataType>(other),
                                                             _extentZero(other._extentZero),
                                                             _extentOne(other._extentOne),
                                                             _extentTwo(other._extentTwo) {}

template <class DataType>
Array3D<DataType>::Array3D(Array3D<DataType> &&other) : _Array_Basic<DataType>(std::move(other)),
                                                        _extentZero(other._extentZero),
                                                        _extentOne(other._extentOne),
                                                        _extentTwo(other._extentTwo) {}

template <class DataType>
Array3D<DataType>::~Array3D() {}

template <class DataType>
Array3D<DataType> &Array3D<DataType>::operator=(Array3D<DataType> const &other)
{
  _Array_Basic<DataType>::operator=(other);
  _extentZero = other._extentZero;
  _extentOne = other._extentOne;
  _extentTwo = other._extentTwo;
  return *this;
}

template <class DataType>
Array3D<DataType> &Array3D<DataType>::operator=(Array3D<DataType> &&other)
{
  _Array_Basic<DataType>::operator=(std::move(other));
  _extentZero = other._extentZero;
  _extentOne = other._extentOne;
  _extentTwo = other._extentTwo;
  return *this;
}

template <class DataType>
inline Array2DView<DataType> const Array3D<DataType>::data_2D(int i) const
{
  return Array2DView<DataType>(_extentOne, _extentTwo, this->m.data() + i * _extentOne * _extentTwo);
}

template <class DataType>
inline Array2DView<DataType> Array3D<DataType>::data_2D(int i)
{
  return Array2DView<DataType>(_extentOne, _extentTwo, this->m.data() + i * _extentOne * _extentTwo);
}

template <class DataType>
inline Array1DView<DataType> const Array3D<DataType>::data_1D(int i, int j) const
{
  return Array1DView<DataType>(_extentTwo, this->m.data() + (i * _extentOne + j) * _extentTwo);
}

template <class DataType>
inline Array1DView<DataType> Array3D<DataType>::data_1D(int i, int j)
{
  return Array1DView<DataType>(_extentTwo, this->m.data() + (i * _extentOne + j) * _extentTwo);
}

template <class DataType>
inline void Array3D<DataType>::resize(int const extentZero,
                                      int const extentOne,
                                      int const extentTwo)
{
  _extentZero = extentZero;
  _extentOne = extentOne;
  _extentTwo = extentTwo;
  std::size_t const _n = _extentZero * _extentOne * _extentTwo;
  this->m.resize(_n);
}

template <class DataType>
inline void Array3D<DataType>::resize(int const extentZero,
                                      int const extentOne,
                                      int const extentTwo,
                                      DataType const new_value)
{
  _extentZero = extentZero;
  _extentOne = extentOne;
  _extentTwo = extentTwo;
  std::size_t const _n = _extentZero * _extentOne * _extentTwo;
  this->m.resize(_n, new_value);
}

template <class DataType>
inline void Array3D<DataType>::resize(int const extentZero,
                                      int const extentOne,
                                      int const extentTwo,
                                      DataType const *new_array)
{
  _extentZero = extentZero;
  _extentOne = extentOne;
  _extentTwo = extentTwo;
  std::size_t const _n = _extentZero * _extentOne * _extentTwo;
  this->m.resize(_n);
  std::copy(new_array, new_array + _n, this->m.data());
}

template <class DataType>
inline const DataType Array3D<DataType>::operator()(int i, int j, int k) const
{
  std::size_t const _n = (i * _extentOne + j) * _extentTwo + k;
  return this->m[_n];
}

template <class DataType>
inline DataType &Array3D<DataType>::operator()(int i, int j, int k)
{
  std::size_t const _n = (i * _extentOne + j) * _extentTwo + k;
  return this->m[_n];
}

template <class DataType>
inline DataType const Array3D<DataType>::at(int i, int j, int k) const
{
  this->_range_check(i, _extentZero);
  this->_range_check(j, _extentOne);
  this->_range_check(k, _extentTwo);
  std::size_t const _n = (i * _extentOne + j) * _extentTwo + k;
  return this->m[_n];
}

template <class DataType>
inline DataType &Array3D<DataType>::at(int i, int j, int k)
{
  this->_range_check(i, _extentZero);
  this->_range_check(j, _extentOne);
  this->_range_check(k, _extentTwo);
  std::size_t const _n = (i * _extentOne + j) * _extentTwo + k;
  return this->m[_n];
}

template <class DataType>
Array3D<DataType>::j_operator::j_operator(Array3D<DataType> &_array, int i) : j_array(_array),
                                                                              _i(i) {}

template <class DataType>
Array3D<DataType>::j_operator::k_operator::k_operator(Array3D<DataType> &_array, int i, int j) : k_array(_array),
                                                                                                 _i(i),
                                                                                                 _j(j) {}

template <class DataType>
DataType const Array3D<DataType>::j_operator::k_operator::operator[](int k) const
{
  std::size_t const _n = (_i * k_array._extentOne + _j) * k_array._extentTwo + k;
  return k_array.m[_n];
}

template <class DataType>
DataType &Array3D<DataType>::j_operator::k_operator::operator[](int k)
{
  std::size_t const _n = (_i * k_array._extentOne + _j) * k_array._extentTwo + k;
  return k_array.m[_n];
}

template <class DataType>
const typename Array3D<DataType>::j_operator::k_operator Array3D<DataType>::j_operator::operator[](int j) const
{
  return k_operator(j_array, _i, j);
}

template <class DataType>
typename Array3D<DataType>::j_operator::k_operator Array3D<DataType>::j_operator::operator[](int j)
{
  return k_operator(j_array, _i, j);
}

template <class DataType>
const typename Array3D<DataType>::j_operator Array3D<DataType>::operator[](int i) const
{
  return j_operator(*this, i);
}

template <class DataType>
typename Array3D<DataType>::j_operator Array3D<DataType>::operator[](int i)
{
  return j_operator(*this, i);
}

#endif // SNAP_HELPER_HPP
