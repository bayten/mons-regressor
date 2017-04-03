#include <algorithm> 

template<class T>
class Vec_t
{
  public:
    T* data;
    int sz;

    Vec_t(int size=0)
    {
        sz = size;
        if(sz > 0)
            data = new T[sz];
        else
            data = NULL;
    }

    Vec_t(int size, T init_vals[])
    {
        sz = size;

        if(sz > 0){
            data = new T[sz];
            for(int i = 0; i < sz; i++)
                data[i] = init_vals[i];
        } else {
            data = NULL;
        }
    }

    ~Vec_t()
    {
        if(data != NULL)
            delete [] data;
    }

    Vec_t(const Vec_t & vec_obj)
    {
        sz = vec_obj.sz;

        data = new T[sz];
        for(int i = 0; i < sz; i++){
            data[i] = vec_obj.data[i];
        }
    }

    Vec_t & operator=(const Vec_t & vec_obj)
    {
        if (this == &vec_obj)  // assigning to myself
            return *this;

        if(data != NULL) {
            delete [] data;
        }

        sz = vec_obj.sz;

        data = new T[sz];
        for(int i = 0; i < sz; i++){
            data[i] = vec_obj.data[i];
        }

        return *this;
    }

    T & operator[](int index)
    {
        if(index == -1){
            return data[sz-1];
        }else if(index > sz){
            std::cout << "ERROR: Index out of bounds!" << std::endl;
            return data[0];
        }
        return data[index];
    }

    void append(T apnd_obj, int index=-1)
    {
        if(data == NULL){
            data = new T[1];
            data[0] = apnd_obj;
            return;
        }

        T* old_data = data;
        data = new T[sz+1];
            
        if(index == -1 || index+1 >= sz){    
            for(int i = 0; i < sz; i++)
                data[i] = old_data[i];

            data[sz] = apnd_obj;

        } else {
            for(int i = 0; i < index; i++)
                data[i] = old_data[i];
            for(int i = index; i < sz; i++)
                data[i+1] = old_data[i+1];
            data[index] = apnd_obj;
        }

        delete [] old_data;
        old_data = NULL;

        sz++;
    }

    Vec_t<T> cutslice(int begin, int end)
    {
        T* old_data = data;
        int slice_sz = end+1-begin;
        data = new T[sz-slice_sz];

        Vec_t<T> carved_slice(slice_sz);

        for(int i = 0; i < begin; i++)
            data[i] = old_data[i];
        for(int i = 0; i < slice_sz; i++)
            carved_slice[i] = old_data[begin+i];
        for(int i = end; i < sz; i++)
            data[i-slice_sz] = old_data[i];

        sz = sz-slice_sz;
        return carved_slice;
    }
};


template<class T>
class Mat_t
{
  public:
    Vec_t< Vec_t<T> > data;
    int sx, sy;

    Mat_t(int x=0, int y=0) : data(x), sx(x), sy(y)
    {
        for(int i = 0; i < sx; i++){
            data[i] = Vec_t<T>(sy);
        }
    }

    Mat_t(int x, int y, T *init_vals[]) : data(x), sx(x), sy(y)
    {
        for(int i = 0; i < sx; i++){
            data[i] = Vec_t<T>(sy);
            for(int j = 0; j < sy; j++)
                data[i][j] = init_vals[i][j];
        }   
    }

    Mat_t(const Mat_t<T> & mat_obj)
    {   
        data = Vec_t< Vec_t<T> >(mat_obj.sx);
        
        sx = mat_obj.sx;
        sy = mat_obj.sy;

        for(int i = 0; i < sx; i++){
            data[i] = mat_obj[i];
        }
    }

    Mat_t(const Vec_t<T> & vec_obj)
    {
        data = Vec_t< Vec_t<T> >(1);

        sx = 1;
        sy = vec_obj.sz;

        data[0] = vec_obj;
    }

    ~Mat_t(){};

    Mat_t<T> & operator=(const Mat_t<T> & mat_obj)
    {
        if (this == &mat_obj)  // assigning to myself
            return *this;

        data = Vec_t< Vec_t<T> >(mat_obj.sx);

        sx = mat_obj.sx;
        sy = mat_obj.sy;

        for(int i = 0; i < sx; i++){
            data[i] = mat_obj[i];
        }

        return *this;
    }

    Vec_t<T> & operator[](int index)
    {
        if(index == -1){
            return data[sx-1];
        }else if(index > sx){
            std::cout << "ERROR: Index out of bounds!" << std::endl;
            return data[0];
        }
        return data[index];
    }
};