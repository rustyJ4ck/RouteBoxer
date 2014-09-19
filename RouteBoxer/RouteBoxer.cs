using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

/**
 * This library is a personal port of RouteBoxer
 * https://code.google.com/p/routeboxer-java/
 * http://google-maps-utility-library-v3.googlecode.com/svn/trunk/routeboxer/src/RouteBoxer.js
 * 
 * @author rustyj4ck 
 * @link https://github.com/rustyJ4ck/RouteBoxer
 */

namespace InforkomApp.RouteBoxer
{
    using RouteBoxerResult = List<RouteBoxer.LatLngBounds>;
    using RouteBoxerInput = List<RouteBoxer.LatLng>;

    public class RouteBoxer
    {
        protected class Logger
        {

            public void warn(params object[] arguments)
            {
                trace(arguments);
            }

            public void trace(params object[] arguments)
            {
                string format = "[RouteBoxer] " + arguments[0].ToString();

                if (arguments.Length > 1)
                {
                    var args = arguments.Skip(1).ToArray();
                    Console.WriteLine(format, args);
                }
                else
                {
                    Console.WriteLine(format);
                }
            }
        }

        protected Logger logger = new Logger();

        private static int EarthRadiusKm = 6371;

        /**
        * @author Cedric NICOLAS
        * see wedrive.mobi for the service using this lib.
        * Utility minimal LatLng class for the need of route boxer 
        */
        public class LatLng //: ICloneable
        {

            public LatLng()
            {
                Lat = 0;
                Lng = 0;

            }

            public LatLng(double lat2, double lng2)
            {
                Lat = lat2;
                Lng = lng2;
            }

            public double Lat
            {
                set;
                get;
            }

            public double Lng
            {
                set;
                get;
            }

            public double LatRad()
            {
                return ToRad(Lat);
            }

            public double LngRad()
            {
                return ToRad(Lng);
            }


            /**
             * A ‘rhumb line’ (or loxodrome) is a path of constant bearing, which crosses all meridians at the same angle.
             * see http://www.movable-type.co.uk/scripts/latlong.html
             * @param brng
             * @param dist
             * @return
             */
            public LatLng RhumbDestinationPoint(double brng, double dist)
            {
                double R = 6378137;
                double d = dist / R;  // d = angular distance covered on earth’s surface
                double lat1 = LatRad(), lon1 = LngRad();
                brng = ToRad(brng);

                double dLat = d * Math.Cos(brng);
                // nasty kludge to overcome ill-conditioned results around parallels of latitude:
                if (Math.Abs(dLat) < 1e-10) dLat = 0; // dLat < 1 mm

                double lat2 = lat1 + dLat;
                double dPhi = Math.Log(Math.Tan(lat2 / 2 + Math.PI / 4) / Math.Tan(lat1 / 2 + Math.PI / 4));
                double q = (dPhi != 0) ? dLat / dPhi : Math.Cos(lat1);  // E-W line gives dPhi=0
                double dLon = d * Math.Sin(brng) / q;

                // check for some daft bugger going past the pole, normalise latitude if so
                if (Math.Abs(lat2) > Math.PI / 2) lat2 = lat2 > 0 ? Math.PI - lat2 : -Math.PI - lat2;

                double lon2 = (lon1 + dLon + 3 * Math.PI) % (2 * Math.PI) - Math.PI;

                return new LatLng(ToDeg(lat2), ToDeg(lon2));
            }

            /**
         * Given a start point and a distance d along constant bearing θ, this will calculate the destination point. 
         * If you maintain a constant bearing along a rhumb line, you will gradually spiral in towards one of the poles.
         * see http://www.movable-type.co.uk/scripts/latlong.html
         * @param dest
         * @return
         */
            public double rhumbBearingTo(LatLng dest)
            {
                double dLon = ToRad(dest.Lng - this.Lng);
                double dPhi = Math.Log(Math.Tan(dest.LatRad() / 2 + Math.PI / 4) / Math.Tan(this.LatRad() / 2 + Math.PI / 4));
                if (Math.Abs(dLon) > Math.PI)
                {
                    dLon = dLon > 0 ? -(2 * Math.PI - dLon) : (2 * Math.PI + dLon);
                }
                return ToBrng(Math.Atan2(dLon, dPhi));
            }

            public override String ToString()
            {
                return FormatLatOrLong(Lat) + "," + FormatLatOrLong(Lng);
            }

            private String FormatLatOrLong(double latOrLng)
            {
                return string.Format("{0:0.00000}", latOrLng);
            }

            public LatLng Clone()
            {
                return new LatLng(Lat, Lng);
            }

            public double DistanceFrom(LatLng otherLatLng)
            {
                double b = Lat * Math.PI / 180.0;
                double c = otherLatLng.Lat * Math.PI / 180.0;
                double d = b - c;
                double e = Lng * Math.PI / 180.0 - otherLatLng.Lng * Math.PI / 180.0;

                double f = 2.0 * Math.Asin(Math.Sqrt(Math.Pow(Math.Sin(d / 2.0), 2.0) + Math.Cos(b) * Math.Cos(c)
                    * Math.Pow(Math.Sin(e / 2.0), 2.0)));
                return f * 6378137.0;
            }
        }

        /**
        * minimal LatLngBounds utility class for RouteBoxer needs
        * i do not guarantee the extend method in southern hemisphere
        */
        public class LatLngBounds
        {
            private LatLng southwest, northeast;

            public LatLngBounds()
            {
            }

            public LatLngBounds(LatLng southwest, LatLng northeast)
            {
                this.southwest = southwest;
                this.northeast = northeast;
            }

            public LatLng SouthWest
            {
                get { return southwest; }
            }

            public LatLng GetSouthWest()
            {
                return southwest;
            }

            public void SetSouthWest(LatLng southwest)
            {
                this.southwest = southwest;
            }

            public LatLng NorthEast
            {
                get { return northeast; }
            }

            public LatLng GetNorthEast()
            {
                return northeast;
            }

            public void SetNorthEast(LatLng northeast)
            {
                this.northeast = northeast;
            }

            public override bool Equals(Object o)
            {
                if (this == o) return true;
                if (o == null || this.GetType().FullName != o.GetType().FullName) return false;

                LatLngBounds that = (LatLngBounds)o;

                if (northeast != null ? !northeast.Equals(that.northeast) : that.northeast != null) return false;
                if (southwest != null ? !southwest.Equals(that.southwest) : that.southwest != null) return false;

                return true;
            }

            public override string ToString()
            {
                return string.Format("{0},{1} - {2},{3}",
                    southwest.Lat, southwest.Lng,
                    northeast.Lat, northeast.Lng
                );
            }

            public void Extend(LatLng latLng)
            {
                if (southwest == null)
                {
                    southwest = latLng.Clone();
                    if (northeast == null) northeast = latLng.Clone();
                    return;
                }
                if (northeast == null)
                {
                    northeast = latLng.Clone();
                    return;
                }

                if (latLng.Lat < southwest.Lat)
                    southwest.Lat = latLng.Lat;
                else if (latLng.Lat > northeast.Lat)
                    northeast.Lat = latLng.Lat;
                if (latLng.Lng < southwest.Lng)
                    southwest.Lng = latLng.Lng;
                else if (latLng.Lng > northeast.Lng)
                    northeast.Lng = latLng.Lng;
            }

            public bool Contains(LatLng latLng)
            {
                if (southwest == null || northeast == null) return false;
                if (latLng.Lat < southwest.Lat) return false;
                else if (latLng.Lat > northeast.Lat) return false;
                else if (latLng.Lng < southwest.Lng) return false;
                else if (latLng.Lng > northeast.Lng) return false;
                return true;
            }

            public LatLng GetCenter()
            {
                return new LatLng(southwest.Lat + (northeast.Lat - southwest.Lat) / 2, southwest.Lng + (northeast.Lng - southwest.Lng) / 2);
            }

            public override int GetHashCode()
            {
                int result = southwest != null ? southwest.GetHashCode() : 0;
                result = 31 * result + (northeast != null ? northeast.GetHashCode() : 0);
                return result;
            }

        }


        private int[][] _grid;
        private List<Double> _latGrid;
        private List<Double> _lngGrid;
        private List<LatLngBounds> _boxesX;
        private List<LatLngBounds> _boxesY;
        /**
         * Creates a new RouteBoxer
         *
         * @constructor
         */
        public RouteBoxer()
        {

        }

        /**
           * utility method to get a LatLng list from an encoded google polyline 
           * You need to pass to box() method below such a long list of LatLngs in order to have it works well
           * @param encodedPoints
           * @return
            */
        public static List<LatLng> DecodePath(String encodedPoints)
        {
            List<LatLng> poly = new List<LatLng>();
            int index = 0, len = encodedPoints.Length;
            int lat = 0, lng = 0;
            while (index < len)
            {
                int b, shift = 0, result = 0;
                do
                {
                    b = encodedPoints.ElementAt(index++) - 63;
                    result |= (b & 0x1f) << shift;
                    shift += 5;
                } while (b >= 0x20);
                int dlat = ((result & 1) != 0 ? ~(result >> 1) : (result >> 1));
                lat += dlat;
                shift = 0;
                result = 0;
                do
                {
                    b = encodedPoints.ElementAt(index++) - 63;
                    result |= (b & 0x1f) << shift;
                    shift += 5;
                } while (b >= 0x20);
                int dlng = ((result & 1) != 0 ? ~(result >> 1) : (result >> 1));
                lng += dlng;
                LatLng p = new LatLng((((double)lat / 1E5)),
                    (((double)lng / 1E5)));
                poly.Add(p);
            }
            return poly;

        }



        /**
         * Generates boxes for a given route and distance
         *
         * @param {LatLng[] } path The path along
         *           which to create boxes. The path object must be either an Array of
         *           LatLng objects
         * @param {Number} range The distance in kms around the route that the generated
         *           boxes must cover.
         * @return {LatLngBounds[]} An List of boxes that covers the whole
         *           path.
         */
        public RouteBoxerResult Boxes(List<LatLng> path, double range)
        {
            // Two dimensional array representing the cells in the grid overlaid on the path
            this._grid = null;

            // Array that holds the latitude coordinate of each vertical grid line
            this._latGrid = new List<Double>();

            // Array that holds the longitude coordinate of each horizontal grid line  
            this._lngGrid = new List<Double>();

            // Array of bounds that cover the whole route formed by merging cells that
            //  the route intersects first horizontally, and then vertically
            this._boxesX = new List<LatLngBounds>();

            // Array of bounds that cover the whole route formed by merging cells that
            //  the route intersects first vertically, and then horizontally
            this._boxesY = new List<LatLngBounds>();

            // The array of LatLngs representing the vertices of the path
            List<LatLng> vertices = null;

            vertices = path;


            // Build the grid that is overlaid on the route
            this._BuildGrid(vertices, range);

            // Identify the grid cells that the route intersects
            this._FindIntersectingCells(vertices);

            // Merge adjacent intersected grid cells (and their neighbours) into two sets
            //  of bounds, both of which cover them completely
            this._MergeIntersectingCells();

            // Return the set of merged bounds that has the fewest elements
            return (this._boxesX.Count <= this._boxesY.Count ?
                this._boxesX :
                this._boxesY);
        }

        private const float GridFactor = 700;

        /**
         * Generates boxes for a given route and distance
         *
         * @param {LatLng[]} vertices The vertices of the path over which to lay the grid
         * @param {Number} range The spacing of the grid cells.
         */
        private void _BuildGrid(List<LatLng> vertices, double range)
        {
            // Create a LatLngBounds object that contains the whole path
            LatLngBounds routeBounds = new LatLngBounds();

            logger.trace("_BuildGrid range: {0}", range);

            range *= GridFactor;

            logger.trace("vertices[0]" + vertices[0].ToString() + " vertices[" + (vertices.Count - 1) + "] " + vertices.ElementAt(vertices.Count - 1).ToString());

            for (int i = 0; i < vertices.Count; i++)
            {
                routeBounds.Extend(vertices[i]);
            }

            logger.trace("routeBounds " + routeBounds.ToString());

            // Find the center of the bounding box of the path
            LatLng routeBoundsCenter = routeBounds.GetCenter();

            logger.trace("routeBoundsCenter " + routeBoundsCenter.ToString());

            // Starting from the center define grid lines outwards vertically until they
            // extend beyond the edge of the bounding box by more than one cell
            this._latGrid.Add(routeBoundsCenter.Lat);
            LatLng rhumb = routeBoundsCenter.RhumbDestinationPoint(0, range);

            logger.trace("rhumb 1 " + rhumb.ToString());

            // Add lines from the center out to the north
            this._latGrid.Add(rhumb.Lat);
            for (int i = 2; this._latGrid[i - 2] < routeBounds.GetNorthEast().Lat; i++)
            {
                this._latGrid.Add(routeBoundsCenter.RhumbDestinationPoint(0, range * i).Lat);
            }

            logger.trace("pass1 latGrid size " + _latGrid.Count);

            // Add lines from the center out to the south  
            for (int i1 = 1; this._latGrid[1] > routeBounds.GetSouthWest().Lat; i1++)
            {
                this._latGrid.Insert(0, routeBoundsCenter.RhumbDestinationPoint(180, range * i1).Lat);
            }

            logger.trace("pass2 latGrid size " + _latGrid.Count);

            // Starting from the center define grid lines outwards horizontally until they
            // extend beyond the edge of the bounding box by more than one cell  
            this._lngGrid.Add(routeBoundsCenter.Lng);

            // Add lines from the center out to the east
            this._lngGrid.Add(routeBoundsCenter.RhumbDestinationPoint(90, range).Lng);
            for (int i2 = 2; this._lngGrid[i2 - 2] < routeBounds.GetNorthEast().Lng; i2++)
            {
                this._lngGrid.Add(routeBoundsCenter.RhumbDestinationPoint(90, range * i2).Lng);
            }

            logger.trace("pass1 lngGrid_ size " + _latGrid.Count);

            // Add lines from the center out to the west
            for (int i3 = 1; this._lngGrid[1] > routeBounds.GetSouthWest().Lng; i3++)
            {
                this._lngGrid.Insert(0, routeBoundsCenter.RhumbDestinationPoint(270, range * i3).Lng);
            }

            logger.trace("pass2 lngGrid_ size " + _latGrid.Count);

            // Create a two dimensional array representing this grid
            // this.grid_ = new int[this.lngGrid_.Count, this.latGrid_.Count];
            this._grid = new int[this._lngGrid.Count][]; // this.latGrid_.Count];

            for (var i = 0; i < this._lngGrid.Count; i++)
            {
                this._grid[i] = new int[this._latGrid.Count];
            }

        }

        /**
         * Find all of the cells in the overlaid grid that the path intersects
         *
         * @param {LatLng[]} vertices The vertices of the path
         */
        private void _FindIntersectingCells(List<LatLng> vertices)
        {
            // Find the cell where the path begins
            int[] hintXY = this.getCellCoords_(vertices[0]);

            // Mark that cell and it's neighbours for inclusion in the boxes
            this._MarkCell(hintXY);

            // Work through each vertex on the path identifying which grid cell it is in
            for (int i = 1; i < vertices.Count; i++)
            {
                // Use the known cell of the previous vertex to help find the cell of this vertex
                int[] gridXY = this.getGridCoordsFromHint_(vertices[i], vertices.ElementAt(i - 1), hintXY);

                if (gridXY[0] == hintXY[0] && gridXY[1] == hintXY[1])
                {
                    // This vertex is in the same cell as the previous vertex
                    // The cell will already have been marked for inclusion in the boxes
                    continue;

                }
                else if ((Math.Abs(hintXY[0] - gridXY[0]) == 1 && hintXY[1] == gridXY[1]) ||
                    (hintXY[0] == gridXY[0] && Math.Abs(hintXY[1] - gridXY[1]) == 1))
                {
                    // This vertex is in a cell that shares an edge with the previous cell
                    // Mark this cell and it's neighbours for inclusion in the boxes
                    this._MarkCell(gridXY);

                }
                else
                {
                    // This vertex is in a cell that does not share an edge with the previous
                    //  cell. This means that the path passes through other cells between
                    //  this vertex and the previous vertex, and we must determine which cells
                    //  it passes through
                    this._GetGridIntersects(vertices[i - 1], vertices[i], hintXY, gridXY);
                }

                // Use this cell to find and compare with the next one
                hintXY = gridXY;
            }
        }

        /**
         * Find the cell a path vertex is in by brute force iteration over the grid
         *
         * @param {LatLng[]} latlng The latlng of the vertex
         * @return {Number[][]} The cell coordinates of this vertex in the grid
         */
        private int[] getCellCoords_(LatLng latlng)
        {
            int x, y;
            for (x = 0; this._lngGrid[x] < latlng.Lng; x++) { }
            for (y = 0; this._latGrid[y] < latlng.Lat; y++) { }
            int[] result = { x - 1, y - 1 };
            return result;
        }

        /**
         * Find the cell a path vertex is in based on the known location of a nearby
         *  vertex. This saves searching the whole grid when working through vertices
         *  on the polyline that are likely to be in close proximity to each other.
         *
         * @param {LatLng[]} latlng The latlng of the vertex to locate in the grid
         * @param {LatLng[]} hintlatlng The latlng of the vertex with a known location
         * @param {Number[]} hint The cell containing the vertex with a known location
         * @return {Number[]} The cell coordinates of the vertex to locate in the grid
         */
        private int[] getGridCoordsFromHint_(LatLng latlng, LatLng hintlatlng, int[] hint)
        {
            int x = 0, y = 0;
            try
            {
                if (latlng.Lng > hintlatlng.Lng)
                {
                    for (x = hint[0]; this._lngGrid[x + 1] < latlng.Lng; x++) { }
                }
                else
                {
                    for (x = hint[0]; this._lngGrid[x] > latlng.Lng; x--) { }
                }

                if (latlng.Lat > hintlatlng.Lat)
                {
                    for (y = hint[1]; this._latGrid[y + 1] < latlng.Lat; y++) { }
                }
                else
                {
                    for (y = hint[1]; this._latGrid[y] > latlng.Lat; y--) { }
                }
            }
            catch (Exception e) //IndexOutOfBoundsException
            {
                logger.warn("getGridCoordsFromHint_ IndexOutOfBoundsException x" + x + " y " + y);
            }
            int[] result = { x, y };
            return result;
        }


        /**
         * Identify the grid squares that a path segment between two vertices
         * intersects with by:
         * 1. Finding the bearing between the start and end of the segment
         * 2. Using the delta between the lat of the start and the lat of each
         *    latGrid boundary to find the distance to each latGrid boundary
         * 3. Finding the lng of the intersection of the line with each latGrid
         *     boundary using the distance to the intersection and bearing of the line
         * 4. Determining the x-coord on the grid of the point of intersection
         * 5. Filling in all squares between the x-coord of the previous intersection
         *     (or start) and the current one (or end) at the current y coordinate,
         *     which is known for the grid line being intersected
         *     
         * @param {LatLng} start The latlng of the vertex at the start of the segment
         * @param {LatLng} end The latlng of the vertex at the end of the segment
         * @param {Number[]} startXY The cell containing the start vertex
         * @param {Number[]} endXY The cell containing the vend vertex
         */
        private void _GetGridIntersects(LatLng start, LatLng end, int[] startXY, int[] endXY)
        {
            LatLng edgePoint;
            int[] edgeXY;
            int i;

            double brng = start.rhumbBearingTo(end);         // Step 1.

            LatLng hint = start;
            int[] hintXY = startXY;

            // Handle a line segment that travels south first
            if (end.Lat > start.Lat)
            {
                // Iterate over the east to west grid lines between the start and end cells
                for (i = startXY[1] + 1; i <= endXY[1]; i++)
                {
                    // Find the latlng of the point where the path segment intersects with
                    //  this grid line (Step 2 & 3)
                    edgePoint = this.getGridIntersect_(start, brng, this._latGrid[i]);

                    // Find the cell containing this intersect point (Step 4)
                    edgeXY = this.getGridCoordsFromHint_(edgePoint, hint, hintXY);

                    // Mark every cell the path has crossed between this grid and the start,
                    //   or the previous east to west grid line it crossed (Step 5)
                    this.fillInGridSquares_(hintXY[0], edgeXY[0], i - 1);

                    // Use the point where it crossed this grid line as the reference for the
                    //  next iteration
                    hint = edgePoint;
                    hintXY = edgeXY;
                }

                // Mark every cell the path has crossed between the last east to west grid
                //  line it crossed and the end (Step 5)
                this.fillInGridSquares_(hintXY[0], endXY[0], i - 1);

            }
            else
            {
                // Iterate over the east to west grid lines between the start and end cells
                for (i = startXY[1]; i > endXY[1]; i--)
                {
                    // Find the latlng of the point where the path segment intersects with
                    //  this grid line (Step 2 & 3)
                    edgePoint = this.getGridIntersect_(start, brng, this._latGrid[i]);

                    // Find the cell containing this intersect point (Step 4)
                    edgeXY = this.getGridCoordsFromHint_(edgePoint, hint, hintXY);

                    // Mark every cell the path has crossed between this grid and the start,
                    //   or the previous east to west grid line it crossed (Step 5)
                    this.fillInGridSquares_(hintXY[0], edgeXY[0], i);

                    // Use the point where it crossed this grid line as the reference for the
                    //  next iteration
                    hint = edgePoint;
                    hintXY = edgeXY;
                }

                // Mark every cell the path has crossed between the last east to west grid
                //  line it crossed and the end (Step 5)
                this.fillInGridSquares_(hintXY[0], endXY[0], i);

            }
        }

        /**
         * Find the latlng at which a path segment intersects with a given
         *   line of latitude
         *     
         * @param {LatLng} start The vertex at the start of the path segment
         * @param {Number} brng The bearing of the line from start to end
         * @param {Number} gridLineLat The latitude of the grid line being intersected
         * @return {LatLng} The latlng of the point where the path segment intersects
         *                    the grid line
         */
        private LatLng getGridIntersect_(LatLng start, double brng, double gridLineLat)
        {
            double d = EarthRadiusKm * ((ToRad(gridLineLat) - start.LatRad()) / Math.Cos(ToRad(brng)));
            return start.RhumbDestinationPoint(brng, d);
        }

        /**
         * Mark all cells in a given row of the grid that lie between two columns
         *   for inclusion in the boxes
         *     
         * @param {Number} startx The first column to include
         * @param {Number} endx The last column to include
         * @param {Number} y The row of the cells to include
         */
        private void fillInGridSquares_(int startx, int endx, int y)
        {
            logger.trace("fillInGridSquares_ startx" + startx + " endx " + endx + " y " + y);
            int x;
            if (startx < endx)
            {
                for (x = startx; x <= endx; x++)
                {
                    int[] cell = { x, y };
                    this._MarkCell(cell);
                }
            }
            else
            {
                for (x = startx; x >= endx; x--)
                {
                    int[] cell = { x, y };
                    this._MarkCell(cell);
                }
            }
        }

        /**
         * Mark a cell and the 8 immediate neighbours for inclusion in the boxes
         *     
         * @param {Number[]} square The cell to mark
         */
        private void _MarkCell(int[] cell)
        {
            int x = cell[0];
            int y = cell[1];
            try
            {
                logger.trace("markCell x" + x + " y " + y);
                this._grid[x - 1][y - 1] = 1;
                this._grid[x][y - 1] = 1;
                this._grid[x + 1][y - 1] = 1;
                this._grid[x - 1][y] = 1;
                this._grid[x][y] = 1;
                this._grid[x + 1][y] = 1;
                this._grid[x - 1][y + 1] = 1;
                this._grid[x][y + 1] = 1;
                this._grid[x + 1][y + 1] = 1;
            }
            catch (Exception e) // IndexOutOfBoundsException
            {
                logger.warn("markCell_ IndexOutOfBoundsException x" + x + " y " + y);
            }
        }

        /**
         * Create two sets of bounding boxes, both of which cover all of the cells that
         *   have been marked for inclusion.
         *
         * The first set is created by combining adjacent cells in the same column into
         *   a set of vertical rectangular boxes, and then combining boxes of the same
         *   height that are adjacent horizontally.
         *
         * The second set is created by combining adjacent cells in the same row into
         *   a set of horizontal rectangular boxes, and then combining boxes of the same
         *   width that are adjacent vertically.
         *     
         */
        private void _MergeIntersectingCells()
        {
            int x, y;
            LatLngBounds box;

            // The box we are currently expanding with new cells
            LatLngBounds currentBox = null;

            // Traverse the grid a row at a time
            for (y = 0; y < this._grid[0].Length; y++)
            {
                for (x = 0; x < this._grid.Length; x++)
                {

                    if (this._grid[x][y] == 1)
                    {
                        // This cell is marked for inclusion. If the previous cell in this
                        //   row was also marked for inclusion, merge this cell into it's box.
                        // Otherwise start a new box.
                        int[] cell = { x, y };
                        box = this._GetCellBounds(cell);
                        if (currentBox != null)
                        {
                            currentBox.Extend(box.GetNorthEast());
                        }
                        else
                        {
                            currentBox = box;
                        }

                    }
                    else
                    {
                        // This cell is not marked for inclusion. If the previous cell was
                        //  marked for inclusion, merge it's box with a box that spans the same
                        //  columns from the row below if possible.
                        this._MergeBoxesY(currentBox);
                        currentBox = null;
                    }
                }
                // If the last cell was marked for inclusion, merge it's box with a matching
                //  box from the row below if possible.
                this._MergeBoxesY(currentBox);
                currentBox = null;
            }

            // Traverse the grid a column at a time
            for (x = 0; x < this._grid.Length; x++)
            {
                for (y = 0; y < this._grid[0].Length; y++)
                {
                    if (this._grid[x][y] == 1)
                    {

                        // This cell is marked for inclusion. If the previous cell in this
                        //   column was also marked for inclusion, merge this cell into it's box.
                        // Otherwise start a new box.
                        int[] cell = { x, y };
                        if (currentBox != null)
                        {

                            box = this._GetCellBounds(cell);
                            currentBox.Extend(box.GetNorthEast());
                        }
                        else
                        {
                            currentBox = this._GetCellBounds(cell);
                        }

                    }
                    else
                    {
                        // This cell is not marked for inclusion. If the previous cell was
                        //  marked for inclusion, merge it's box with a box that spans the same
                        //  rows from the column to the left if possible.
                        this._MergeBoxesX(currentBox);
                        currentBox = null;

                    }
                }
                // If the last cell was marked for inclusion, merge it's box with a matching
                //  box from the column to the left if possible.
                this._MergeBoxesX(currentBox);
                currentBox = null;
            }
        }

        /**
         * Search for an existing box in an adjacent row to the given box that spans the
         * same set of columns and if one is found merge the given box into it. If one
         * is not found, append this box to the list of existing boxes.
         *
         * @param {LatLngBounds}  The box to merge
         */
        private void _MergeBoxesX(LatLngBounds box)
        {
            if (box != null)
            {
                for (int i = 0; i < this._boxesX.Count; i++)
                {
                    if (Math.Abs(this._boxesX[i].GetNorthEast().Lng - box.GetSouthWest().Lng) < 0.001 &&
                        Math.Abs(this._boxesX[i].GetSouthWest().Lat - box.GetSouthWest().Lat) < 0.001 &&
                        Math.Abs(this._boxesX[i].GetNorthEast().Lat - box.GetNorthEast().Lat) < 0.001)
                    {
                        this._boxesX[i].Extend(box.GetNorthEast());
                        return;
                    }
                }
                this._boxesX.Add(box);
            }
        }

        /**
         * Search for an existing box in an adjacent column to the given box that spans
         * the same set of rows and if one is found merge the given box into it. If one
         * is not found, append this box to the list of existing boxes.
         *
         * @param {LatLngBounds}  The box to merge
         */
        private void _MergeBoxesY(LatLngBounds box)
        {
            if (box != null)
            {
                for (int i = 0; i < this._boxesY.Count; i++)
                {
                    if (Math.Abs(this._boxesY[i].GetNorthEast().Lat - box.GetSouthWest().Lat) < 0.001 &&
                        Math.Abs(this._boxesY[i].GetSouthWest().Lng - box.GetSouthWest().Lng) < 0.001 &&
                        Math.Abs(this._boxesY[i].GetNorthEast().Lng - box.GetNorthEast().Lng) < 0.001)
                    {
                        this._boxesY[i].Extend(box.GetNorthEast());
                        return;
                    }
                }
                this._boxesY.Add(box);
            }
        }

        /**
         * Obtain the LatLng of the origin of a cell on the grid
         *
         * @param {Number[]} cell The cell to lookup.
         * @return {LatLng} The latlng of the origin of the cell.
         */
        private LatLngBounds _GetCellBounds(int[] cell)
        {
            return new LatLngBounds(
                new LatLng(this._latGrid.ElementAt(cell[1]), this._lngGrid.ElementAt(cell[0])),
                new LatLng(this._latGrid.ElementAt(cell[1] + 1), this._lngGrid.ElementAt(cell[0] + 1)));
        }



        /**
         * Extend the Number object to convert degrees to radians
         *
         * @return {Number} Bearing in radians
         * @ignore
         */
        public static double ToRad(double value)
        {
            return value * Math.PI / 180;
        }

        /**
         * Extend the Number object to convert radians to degrees
         *
         * @return {Number} Bearing in degrees
         * @ignore
         */
        public static double ToDeg(double value)
        {
            return value * 180 / Math.PI;
        }

        /**
         * Normalize a heading in degrees to between 0 and +360
         *
         * @return {Number} Return 
         * @ignore
         */
        public static double ToBrng(double value)
        {
            return (ToDeg(value) + 360) % 360;
        }
    }


}

