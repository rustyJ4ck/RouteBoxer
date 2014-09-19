private void DrawRouteBoxes(List<RouteBoxer.RouteBoxer.LatLng> steps) 
        {

            var boxes = new RouteBoxer.RouteBoxer().Boxes(steps, 5);

            boxes.ForEach(box => {

                var po = new PolygonOptions();

                po.Add(
                    new LatLng(box.SouthWest.Lat, box.SouthWest.Lng), 
                    new LatLng(box.SouthWest.Lat, box.NorthEast.Lng), 
                    new LatLng(box.NorthEast.Lat, box.NorthEast.Lng), 
                    new LatLng(box.NorthEast.Lat, box.SouthWest.Lng), 
                    new LatLng(box.SouthWest.Lat, box.SouthWest.Lng)
                );

                po.InvokeStrokeColor(Color.PaleVioletRed);
                po.InvokeFillColor(Color.Transparent);
                po.InvokeStrokeWidth(2.0f);

                _map.AddPolygon(po);

            });
        }