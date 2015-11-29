/****************************************************************************
**
** For Copyright & Licensing information, see COPYRIGHT in project root
**
****************************************************************************/

#include <QImage>
#include <QVector>

typedef QVector<unsigned int> IntArray;
typedef QVector<QPoint>       PointArray;
typedef QVector<IntArray>     Image;

class HoughCircleDetector
{
  public: /* class */
  
    HoughCircleDetector() {}
   ~HoughCircleDetector() {}
  
  public: /* methods */
  
    QImage detect(const QImage &source, unsigned int min_r, unsigned int max_r);
  
  private: /* methods */
    void accum_circle(Image &image, const QPoint &position, unsigned int radius);
    void accum_circle2(Image &image, const QPoint &p1, const QPoint &p2, unsigned int radius);
    void accum_circle4(Image &image, const QPoint &p1, const QPoint &p2, const QPoint &p3, const QPoint &p4, unsigned int radius);

    void accum_circle_row(Image &image, unsigned int row, const IntArray &col_indices, unsigned int radius);

    void draw_circle(QImage &image, const QPoint &position, unsigned int radius, const QColor &color);
    void draw_pixel(QImage &image, const int x, const int y, const QRgb &rgb);
    
    QImage edges(const QImage &source);
};

