using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using InforkomApp.RouteBoxer;

namespace TestRouteBoxer
{
    class Program
    {
        static void Main(string[] args)
        {
            
            string points = @"qdjqIyupcFDq@NmDh@sKJeCDy@B_@L_CXkGFgB@_@@QAW?A?QASCSCSCUEMAEESGOKYIQIMKMMQm@m@IIoAgAs@i@UOOMa@Og@Q]McBe@MCgBg@aAWmA]kA_@m@YWQKMIMQYSiA@qCAyCDuB?s@Ea@GYIWO[_@{@c@}@i@oA}@oBmBgESk@Ia@S{@g@}CMcAKq@]uCOq@G[E_@E]IuBEsACc@Gc@E_@Gc@EYCO?Q?o@FwAHkCFgBL{D@SFoAB]H{ABYFwALaCNqCReD";

            var decodedPoints = RouteBoxer.DecodePath(points);

            var boxer = new InforkomApp.RouteBoxer.RouteBoxer();

            var results = boxer.Boxes(decodedPoints, 10);

            results.ForEach(s =>
            {
                Console.WriteLine(s.ToString());
            });
            

            Console.ReadKey();
        }
    }
}
