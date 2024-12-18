/** \file       hoMRImage.hxx
    \brief      Implementation of N-dimensional MR image class for gadgetron
    \author     Hui Xue
*/

namespace Gadgetron
{
    template <typename T, unsigned int D>
    hoMRImage<T, D>::hoMRImage () : BaseClass()
    {
        memset(&this->header_, 0, sizeof(mrd::ImageHeader));
    }

    template <typename T, unsigned int D>
    hoMRImage<T, D>::hoMRImage (const std::vector<size_t>& dimensions) : BaseClass( const_cast<std::vector<size_t>& >(dimensions) )
    {
        this->create(dimensions);
    }

    template <typename T, unsigned int D>
    hoMRImage<T, D>::hoMRImage (const std::vector<size_t>& dimensions,
        const std::vector<coord_type>& pixelSize) : BaseClass( const_cast<std::vector<size_t>& >(dimensions) )
    {
        this->create(dimensions, pixelSize);
    }

    template <typename T, unsigned int D>
    hoMRImage<T, D>::hoMRImage (const std::vector<size_t>& dimensions,
        const std::vector<coord_type>& pixelSize,
        const std::vector<coord_type>& origin) : BaseClass( const_cast<std::vector<size_t>& >(dimensions) )
    {
        this->create(dimensions, pixelSize, origin);
    }

    template <typename T, unsigned int D>
    hoMRImage<T, D>::hoMRImage (const std::vector<size_t>& dimensions,
                                const std::vector<coord_type>& pixelSize,
                                const std::vector<coord_type>& origin,
                                const axis_type& axis) : BaseClass( const_cast<std::vector<size_t>& >(dimensions) )
    {
        this->create(dimensions, pixelSize, origin, axis);
    }

    template <typename T, unsigned int D>
    hoMRImage<T, D>::hoMRImage(size_t len) : BaseClass(len)
    {
        std::vector<size_t> dimension(1, len);
        this->create(dimension);
    }

    template <typename T, unsigned int D>
    hoMRImage<T, D>::hoMRImage(size_t sx, size_t sy) : BaseClass(sx, sy)
    {
        std::vector<size_t> dimension(2);
        dimension[0] = sx;
        dimension[1] = sy;
        this->create(dimension);
    }

    template <typename T, unsigned int D>
    hoMRImage<T, D>::hoMRImage(size_t sx, size_t sy, size_t sz) : BaseClass(sx, sy, sz)
    {
        std::vector<size_t> dimension(3);
        dimension[0] = sx;
        dimension[1] = sy;
        dimension[2] = sz;
        this->create(dimension);
    }

    template <typename T, unsigned int D>
    hoMRImage<T, D>::hoMRImage(size_t sx, size_t sy, size_t sz, size_t st) : BaseClass(sx, sy, sz, st)
    {
        std::vector<size_t> dimension(4);
        dimension[0] = sx;
        dimension[1] = sy;
        dimension[2] = sz;
        dimension[3] = st;
        this->create(dimension);
    }

    template <typename T, unsigned int D>
    hoMRImage<T, D>::hoMRImage(size_t sx, size_t sy, size_t sz, size_t st, size_t sp) : BaseClass(sx, sy, sz, st, sp)
    {
        std::vector<size_t> dimension(5);
        dimension[0] = sx;
        dimension[1] = sy;
        dimension[2] = sz;
        dimension[3] = st;
        dimension[4] = sp;
        this->create(dimension);
    }

    template <typename T, unsigned int D>
    hoMRImage<T, D>::hoMRImage(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq) : BaseClass(sx, sy, sz, st, sp, sq)
    {
        std::vector<size_t> dimension(6);
        dimension[0] = sx;
        dimension[1] = sy;
        dimension[2] = sz;
        dimension[3] = st;
        dimension[4] = sp;
        dimension[5] = sq;
        this->create(dimension);
    }

    template <typename T, unsigned int D>
    hoMRImage<T, D>::hoMRImage(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, size_t sr) : BaseClass(sx, sy, sz, st, sp, sq, sr)
    {
        std::vector<size_t> dimension(7);
        dimension[0] = sx;
        dimension[1] = sy;
        dimension[2] = sz;
        dimension[3] = st;
        dimension[4] = sp;
        dimension[5] = sq;
        dimension[6] = sr;
        this->create(dimension);
    }

    template <typename T, unsigned int D>
    hoMRImage<T, D>::hoMRImage(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, size_t sr, size_t ss) : BaseClass(sx, sy, sz, st, sp, sq, sr, ss)
    {
        std::vector<size_t> dimension(8);
        dimension[0] = sx;
        dimension[1] = sy;
        dimension[2] = sz;
        dimension[3] = st;
        dimension[4] = sp;
        dimension[5] = sq;
        dimension[6] = sr;
        dimension[7] = ss;
        this->create(dimension);
    }

    template <typename T, unsigned int D>
    hoMRImage<T, D>::hoMRImage (const std::vector<size_t>& dimensions, T* data, bool delete_data_on_destruct) : BaseClass(dimensions, data, delete_data_on_destruct)
    {
        this->create(dimensions, data, delete_data_on_destruct);
    }

    template <typename T, unsigned int D>
    hoMRImage<T, D>::hoMRImage (const std::vector<size_t>& dimensions, const std::vector<coord_type>& pixelSize, T* data, bool delete_data_on_destruct) : BaseClass(dimensions, data, delete_data_on_destruct)
    {
        this->create(dimensions, pixelSize, data, delete_data_on_destruct);
    }

    template <typename T, unsigned int D>
    hoMRImage<T, D>::hoMRImage (const std::vector<size_t>& dimensions, const std::vector<coord_type>& pixelSize, const std::vector<coord_type>& origin, T* data, bool delete_data_on_destruct) : BaseClass(dimensions, data, delete_data_on_destruct)
    {
        this->create(dimensions, pixelSize, origin, data, delete_data_on_destruct);
    }

    template <typename T, unsigned int D>
    hoMRImage<T, D>::hoMRImage (const std::vector<size_t>& dimensions, const std::vector<coord_type>& pixelSize, const std::vector<coord_type>& origin, const axis_type& axis, T* data, bool delete_data_on_destruct) : BaseClass(dimensions, data, delete_data_on_destruct)
    {
        this->create(dimensions, pixelSize, origin, axis, data, delete_data_on_destruct);
    }

    template <typename T, unsigned int D>
    hoMRImage<T, D>::hoMRImage(size_t len, T* data, bool delete_data_on_destruct) : BaseClass(len, data, delete_data_on_destruct)
    {
        std::vector<size_t> dimension(1, len);
        this->create(dimension, data, delete_data_on_destruct);
    }

    template <typename T, unsigned int D>
    hoMRImage<T, D>::hoMRImage(size_t sx, size_t sy, T* data, bool delete_data_on_destruct) : BaseClass(sx, sy, data, delete_data_on_destruct)
    {
        std::vector<size_t> dimension(2);
        dimension[0] = sx;
        dimension[1] = sy;
        this->create(dimension, data, delete_data_on_destruct);
    }

    template <typename T, unsigned int D>
    hoMRImage<T, D>::hoMRImage(size_t sx, size_t sy, size_t sz, T* data, bool delete_data_on_destruct) : BaseClass(sx, sy, sz, data, delete_data_on_destruct)
    {
        std::vector<size_t> dimension(3);
        dimension[0] = sx;
        dimension[1] = sy;
        dimension[2] = sz;
        this->create(dimension, data, delete_data_on_destruct);
    }

    template <typename T, unsigned int D>
    hoMRImage<T, D>::hoMRImage(size_t sx, size_t sy, size_t sz, size_t st, T* data, bool delete_data_on_destruct) : BaseClass(sx, sy, sz, st, data, delete_data_on_destruct)
    {
        std::vector<size_t> dimension(4);
        dimension[0] = sx;
        dimension[1] = sy;
        dimension[2] = sz;
        dimension[3] = st;
        this->create(dimension, data, delete_data_on_destruct);
    }

    template <typename T, unsigned int D>
    hoMRImage<T, D>::hoMRImage(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, T* data, bool delete_data_on_destruct) : BaseClass(sx, sy, sz, st, sp, data, delete_data_on_destruct)
    {
        std::vector<size_t> dimension(5);
        dimension[0] = sx;
        dimension[1] = sy;
        dimension[2] = sz;
        dimension[3] = st;
        dimension[4] = sp;
        this->create(dimension, data, delete_data_on_destruct);
    }

    template <typename T, unsigned int D>
    hoMRImage<T, D>::hoMRImage(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, T* data, bool delete_data_on_destruct) : BaseClass(sx, sy, sz, st, sp, sq, data, delete_data_on_destruct)
    {
        std::vector<size_t> dimension(6);
        dimension[0] = sx;
        dimension[1] = sy;
        dimension[2] = sz;
        dimension[3] = st;
        dimension[4] = sp;
        dimension[5] = sq;
        this->create(dimension, data, delete_data_on_destruct);
    }

    template <typename T, unsigned int D>
    hoMRImage<T, D>::hoMRImage(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, size_t sr, T* data, bool delete_data_on_destruct) : BaseClass(sx, sy, sz, st, sp, sq, sr, data, delete_data_on_destruct)
    {
        std::vector<size_t> dimension(7);
        dimension[0] = sx;
        dimension[1] = sy;
        dimension[2] = sz;
        dimension[3] = st;
        dimension[4] = sp;
        dimension[5] = sq;
        dimension[6] = sr;
        this->create(dimension, data, delete_data_on_destruct);
    }

    template <typename T, unsigned int D>
    hoMRImage<T, D>::hoMRImage(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, size_t sr, size_t ss, T* data, bool delete_data_on_destruct) : BaseClass(sx, sy, sz, st, sp, sq, sr, ss, data, delete_data_on_destruct)
    {
        std::vector<size_t> dimension(8);
        dimension[0] = sx;
        dimension[1] = sy;
        dimension[2] = sz;
        dimension[3] = st;
        dimension[4] = sp;
        dimension[5] = sq;
        dimension[6] = sr;
        dimension[7] = ss;
        this->create(dimension, data, delete_data_on_destruct);
    }

    template <typename T, unsigned int D>
    hoMRImage<T, D>::hoMRImage(const hoNDArray<T>& a) : BaseClass(a)
    {
         this->create(a.get_dimensions());
         memcpy(this->data_, a.begin(), this->get_number_of_bytes());
    }

    template <typename T, unsigned int D>
    hoMRImage<T, D>::hoMRImage(const Self& a) : BaseClass()
    {
        *this = a;
    }

    template <typename T, unsigned int D>
    hoMRImage<T, D>& hoMRImage<T, D>::operator=(const Self& rhs)
    {
        if ( &rhs == this ) return *this;

        BaseClass::operator=(rhs);

        this->header_ = rhs.header_;
        this->attrib_ = rhs.attrib_;

        return *this;
    }

    template <typename T, unsigned int D>
    hoMRImage<T, D>::~hoMRImage()
    {
        if (this->delete_data_on_destruct_)
        {
            this->deallocate_memory();
        }
    }

    template <typename T, unsigned int D>
    void hoMRImage<T, D>::clear()
    {
        BaseClass::clear();
        memset(&this->header_, 0, sizeof(mrd::ImageHeader));
        this->attrib_ = mrd::ImageMeta();
    }

    template <typename T, unsigned int D>
    void hoMRImage<T, D>::create(const std::vector<size_t>& dimensions)
    {
        BaseClass::create(dimensions);
        memset(&this->header_, 0, sizeof(mrd::ImageHeader));
    }

    template <typename T, unsigned int D>
    void hoMRImage<T, D>::create(const std::vector<size_t>& dimensions, const std::vector<coord_type>& pixelSize)
    {
        BaseClass::create(dimensions, pixelSize);
        memset(&this->header_, 0, sizeof(mrd::ImageHeader));
    }

    template <typename T, unsigned int D>
    void hoMRImage<T, D>::create(const std::vector<size_t>& dimensions, const std::vector<coord_type>& pixelSize, const std::vector<coord_type>& origin)
    {
        BaseClass::create(dimensions, pixelSize, origin);
        memset(&this->header_, 0, sizeof(mrd::ImageHeader));
    }

    template <typename T, unsigned int D>
    void hoMRImage<T, D>::create(const std::vector<size_t>& dimensions, const std::vector<coord_type>& pixelSize, const std::vector<coord_type>& origin, const axis_type& axis)
    {
        BaseClass::create(dimensions, pixelSize, origin, axis);
        memset(&this->header_, 0, sizeof(mrd::ImageHeader));
    }

    template <typename T, unsigned int D>
    void hoMRImage<T, D>::create(const std::vector<size_t>& dimensions, T* data, bool delete_data_on_destruct)
    {
        BaseClass::create(dimensions, data, delete_data_on_destruct);
        memset(&this->header_, 0, sizeof(mrd::ImageHeader));
    }

    template <typename T, unsigned int D>
    void hoMRImage<T, D>::create(const std::vector<size_t>& dimensions, const std::vector<coord_type>& pixelSize, T* data, bool delete_data_on_destruct)
    {
        BaseClass::create(dimensions, pixelSize, data, delete_data_on_destruct);
        memset(&this->header_, 0, sizeof(mrd::ImageHeader));
    }

    template <typename T, unsigned int D>
    void hoMRImage<T, D>::create(const std::vector<size_t>& dimensions, const std::vector<coord_type>& pixelSize, const std::vector<coord_type>& origin, T* data, bool delete_data_on_destruct)
    {
        BaseClass::create(dimensions, pixelSize, origin, data, delete_data_on_destruct);
        memset(&this->header_, 0, sizeof(mrd::ImageHeader));
    }

    template <typename T, unsigned int D>
    void hoMRImage<T, D>::create(const std::vector<size_t>& dimensions, const std::vector<coord_type>& pixelSize, const std::vector<coord_type>& origin, const axis_type& axis, T* data, bool delete_data_on_destruct)
    {
        BaseClass::create(dimensions, pixelSize, origin, axis, data, delete_data_on_destruct);
        memset(&this->header_, 0, sizeof(mrd::ImageHeader));
    }

    template <typename T, unsigned int D>
    void hoMRImage<T, D>::get_sub_image(const std::vector<size_t>& start, std::vector<size_t>& size, Self& out)
    {
        BaseClass::get_sub_image(start, size, out);

        out.header_ = this->header_;
        out.attrib_ = this->attrib_;
    }

    template <typename T, unsigned int D>
    inline void hoMRImage<T, D>::printContent(std::ostream& os) const
    {
        using namespace std;
        BaseClass::printContent(os);
    }
}
