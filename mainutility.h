#pragma once

template <class T>
class MyVector : public std::vector<T>
{
private:
	int numDims;

public:

	MyVector(int dim)
	{
		this->numDims = dim;
	}

	int getDims()
	{
		return this->numDims;
	}

	inline size_t kdtree_get_point_count() const
	{
		return this->size() / numDims;
	}

	inline T kdtree_get_pt(const size_t idx, const size_t dim) const 
	{
		int index = numDims * idx + dim;
		assert(index < this->size());

		return this->operator[](index);
	}

	template <class BBOX> 
	bool kdtree_get_bbox(BBOX& bb) const
	{
		return false;
	}
};

struct Quadrante
{
	double xMin;
	double xMax;
	double yMin;
	double yMax;

	bool isInside(double x, double y)
	{
		return (xMin <= x && x <= xMax) && (yMin <= y && y <= yMax);
	}
};

struct XYZ
{
	double x, y, z;
};

typedef struct {
	double r;       // a fraction between 0 and 1
	double g;       // a fraction between 0 and 1
	double b;       // a fraction between 0 and 1
} rgb;

typedef struct {
	double h;       // angle in degrees
	double s;       // a fraction between 0 and 1
	double v;       // a fraction between 0 and 1
} hsv;

//da tinta a rgb
rgb hsv2rgb(hsv in)
{
	double      hh, p, q, t, ff;
	long        i;
	rgb         out;

	if (in.s <= 0.0) {       // < is bogus, just shuts up warnings
		out.r = in.v;
		out.g = in.v;
		out.b = in.v;
		return out;
	}
	hh = in.h;
	if (hh >= 360.0) hh = 0.0;
	hh /= 60.0;
	i = (long)hh;
	ff = hh - i;
	p = in.v * (1.0 - in.s);
	q = in.v * (1.0 - (in.s * ff));
	t = in.v * (1.0 - (in.s * (1.0 - ff)));

	switch (i) {
	case 0:
		out.r = in.v;
		out.g = t;
		out.b = p;
		break;
	case 1:
		out.r = q;
		out.g = in.v;
		out.b = p;
		break;
	case 2:
		out.r = p;
		out.g = in.v;
		out.b = t;
		break;

	case 3:
		out.r = p;
		out.g = q;
		out.b = in.v;
		break;
	case 4:
		out.r = t;
		out.g = p;
		out.b = in.v;
		break;
	case 5:
	default:
		out.r = in.v;
		out.g = p;
		out.b = q;
		break;
	}
	return out;
}


struct PaddingVector 
{
	char padding[64]; //padding per evitare false sharing
	std::vector<double> arr;
};
