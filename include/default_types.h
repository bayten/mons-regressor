/* Copyright 2017 Baytekov Nikita */
#ifndef INCLUDE_DEFAULT_TYPES_H_
#define INCLUDE_DEFAULT_TYPES_H_

#include <iostream>
#include <algorithm>
#include <vector>

template<typename T>
class Vec {
 private:
    T* data;
    int sz;
    int real_sz;

 public:
    explicit Vec(int size = 0);  // default constructor
    Vec(int size, T init_vals[]);  // -//-
    Vec(const Vec<T>& vec_obj);  // copy constructor
    explicit Vec(Vec<T>&& vec_obj);  // move constructor
    ~Vec();  // destructor

    Vec<T>& operator=(const Vec<T>& vec_obj);  // copy assignment
    Vec<T>& operator=(Vec<T>&& vec_obj);  // move assignment
    T& operator[](int index);
    const T& operator[](int index) const;

    void append(const T& apnd_obj, int index = -1);
    bool erase(int index = -1);
    Vec<T> cutslice(int begin, int end) const;
    Vec<T> sort() const;
    Vec<int> sort_indices() const;
    int get_size() const {
        return sz;
    }
};


template<typename T>
class Mat {
 private:
    Vec< Vec<T> > data;
    int sx, sy;

 public:
    explicit Mat(int x = 0, int y = 0);
    Mat(int x, int y, T *init_vals[]);
    Mat(const Mat<T>& mat_obj);
    explicit Mat(const Vec<T> & vec_obj);
    ~Mat() {}

    Mat<T>& operator=(const Mat<T>& mat_obj);
    Vec<T>& operator[](int index);
    const Vec<T>& operator[](int index) const;

    int get_sx() const { return sx; }
    int get_sy() const { return sy; }
};


template<typename T>
Vec<T>::Vec(int size) {
    sz = size;
    real_sz = size;

    if (sz > 0)
        data = new T[sz];
    else
        data = nullptr;
}

template<typename T>
Vec<T>::Vec(int size, T init_vals[]) {
    sz = size;
    real_sz = size;

    if (sz > 0) {
        data = new T[sz];
        for (int i = 0; i < sz; i++)
            data[i] = init_vals[i];
    } else {
        data = nullptr;
    }
}

template<typename T>
Vec<T>::Vec(const Vec<T>& vec_obj) {
    sz = vec_obj.sz;
    real_sz = vec_obj.real_sz;
    data = new T[real_sz];
    for (int i = 0; i < sz; i++)
        data[i] = vec_obj.data[i];
}

template<typename T>
Vec<T>::Vec(Vec<T>&& vec_obj) {
    sz = vec_obj.sz;
    real_sz = vec_obj.real_sz;
    data = vec_obj.data;
    vec_obj.data = nullptr;
}

template<typename T>
Vec<T>::~Vec() {
    if (data != nullptr)
        delete [] data;
}

template<typename T>
Vec<T>& Vec<T>::operator=(const Vec<T>& vec_obj) {
    if (this == &vec_obj)  // assigning to myself
        return *this;

    if (data != nullptr)
        delete [] data;

    sz = vec_obj.sz;
    real_sz = vec_obj.real_sz;

    data = new T[real_sz];
    for (int i = 0; i < sz; i++)
        data[i] = vec_obj.data[i];

    return *this;
}

template<typename T>
Vec<T>& Vec<T>::operator=(Vec<T>&& vec_obj) {
    T* tmp_data = data;
    sz = vec_obj.sz;
    real_sz = vec_obj.real_sz;
    data = vec_obj.data;

    vec_obj.data = tmp_data;
    tmp_data = nullptr;

    return (*this);
}

template<typename T>
T& Vec<T>::operator[](int index) {
    if (index == -1) {
        return data[sz-1];
    } else if (index > sz) {
        std::cout << "ERROR: Index out of bounds!" << std::endl;
        return data[0];
    }
    return data[index];
}

template<typename T>
const T& Vec<T>::operator[](int index) const {
    if (index == -1) {
        return data[sz-1];
    } else if (index > sz) {
        std::cout << "ERROR: Index out of bounds!" << std::endl;
        return data[0];
    }
    return data[index];
}

template<typename T>
void Vec<T>::append(const T& apnd_obj, int index) {
    if (data == nullptr) {
        data = new T[1];
        data[0] = apnd_obj;
        sz = 1;
        real_sz = 1;
        return;
    } else if (sz+1 <= real_sz) {
        for (int i = sz; i > index; i--)
            data[i] = data[i-1];
        data[index] = apnd_obj;
        sz++;
        return;
    } else {
        real_sz *= 2;
        T* old_data = data;
        data = new T[real_sz];

        for (int i = 0; i < index; i++)
            data[i] = old_data[i];
        for (int i = index; i < sz; i++)
            data[i+1] = old_data[i+1];
        data[index] = apnd_obj;

        delete [] old_data;
        old_data = nullptr;

        sz++;
    }
}

template<typename T>
bool Vec<T>::erase(int index) {
    if (data == nullptr || (index < 0 && index != -1) || index >= sz)
        return 0;

    for (int i = index; i < sz-1; i++)
            data[i] = data[i+1];
    sz--;

    return 1;
}

template<typename T>
Vec<T> Vec<T>::cutslice(int begin, int end) const {
    int slice_sz = end+1-begin;
    Vec<T> carved_slice(slice_sz);

    for (int i = 0; i < slice_sz; i++)
        carved_slice[i] = data[begin+i];
    for (int i = end; i < sz; i++)
        data[i-slice_sz] = data[i];

    sz = sz-slice_sz;
    return carved_slice;
}

template<typename T>
Vec<T> Vec<T>::sort() const {
    T* new_data = new T[sz];
    for (int i = 0; i < sz; i++)
        new_data[i] = data[i];

    std::sort(new_data, new_data+sz);
    return Vec<T>(sz, new_data);
}

template<typename T>
Vec<int> Vec<T>::sort_indices() const {
    std::vector<int>idx(sz);
    std::iota(std::begin(idx), std::end(idx), 0);
    std::sort(std::begin(idx), std::end(idx),
              [this](int i1, int i2) { return this->data[i1] < this->data[i2]; });
    return Vec<int>(sz, idx);
}


template<typename T>
Mat<T>::Mat(int x, int y) : data(x), sx(x), sy(y) {
    for (int i = 0; i < sx; i++) {
        data[i] = Vec<T>(sy);
    }
}

template<typename T>
Mat<T>::Mat(int x, int y, T *init_vals[]) : data(x), sx(x), sy(y) {
    for (int i = 0; i < sx; i++) {
        data[i] = Vec<T>(sy);
        for (int j = 0; j < sy; j++)
            data[i][j] = init_vals[i][j];
    }
}

template<typename T>
Mat<T>::Mat(const Mat<T>& mat_obj) : data(Vec< Vec<T> >(mat_obj.sx)) {
    sx = mat_obj.sx;
    sy = mat_obj.sy;

    for (int i = 0; i < sx; i++) {
        data[i] = mat_obj[i];
    }
}

template<typename T>
Mat<T>::Mat(const Vec<T> & vec_obj) : data(Vec< Vec<T> >(1)) {
    sx = 1;
    sy = vec_obj.sz;

    data[0] = vec_obj;
}

template<typename T>
Mat<T>& Mat<T>::operator=(const Mat<T> & mat_obj) {
    if (this == &mat_obj)  // assigning to myself
        return *this;

    data = Vec< Vec<T> >(mat_obj.sx);

    sx = mat_obj.sx;
    sy = mat_obj.sy;

    for (int i = 0; i < sx; i++) {
        data[i] = mat_obj[i];
    }

    return *this;
}

template<typename T>
Vec<T>& Mat<T>::operator[](int index) {
    if (index == -1) {
        return data[sx-1];
    } else if (index > sx) {
        std::cout << "ERROR: Index out of bounds!" << std::endl;
        return data[0];
    }
    return data[index];
}

template<typename T>
const Vec<T>& Mat<T>::operator[](int index) const {
    if (index == -1) {
        return data[sx-1];
    } else if (index > sx) {
        std::cout << "ERROR: Index out of bounds!" << std::endl;
        return data[0];
    }
    return data[index];
}

#endif  // INCLUDE_DEFAULT_TYPES_H_
